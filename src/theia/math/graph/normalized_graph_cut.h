// Copyright (C) 2015 The Regents of the University of California (Regents).
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//
//     * Neither the name of The Regents or University of California nor the
//       names of its contributors may be used to endorse or promote products
//       derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// Please contact the author of this library if you have any questions.
// Author: Chris Sweeney (cmsweeney@cs.ucsb.edu)

#ifndef THEIA_MATH_GRAPH_NORMALIZED_GRAPH_CUT_H_
#define THEIA_MATH_GRAPH_NORMALIZED_GRAPH_CUT_H_

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <glog/logging.h>

#include <algorithm>
#include <limits>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "spectra/include/SymEigsSolver.h"

#include "theia/math/matrix/linear_operator.h"
#include "theia/util/hash.h"
#include "theia/util/map_util.h"

namespace theia {

// Given a weighted graph G = {V, E} and edge weights W, the normalized graph
// cut algorithm computes a cut through the graph to segment the graph into two
// subgraphs based on their connectivity and edge weights. A simple, efficient
// algorithm is proposed in: "Normalized Cuts and Image Segmentation" by Shi and
// Malik (PAMI 2000) that solves for the cut by sparse eigen-decomposition.
//
// The input is a set of undirected edges with nodes of type T and double-type
// weights. The two subgraphs are returned with each vector containing node
// identifiers. Each node participates in exactly one of the
// sub-graphs. Additionally, the cost of the cut is an optional output (pass in
// NULL if the cost is not desired) and could be used to determine the stability
// of the cut.
template <typename T>
class NormalizedGraphCut {
 public:
  struct Options {
    // After doing some math, the algorithm will test several points to try and
    // make a cut. The point that has the lowest normalized cut cost will then
    // be used to make the cut. This parameter controls how many points to test
    // when determining the cutting point.
    int num_cuts_to_test = 20;
  };

  explicit NormalizedGraphCut(const Options& options) : options_(options) {}

  // Computes a graph cut and optionally returns the cost of the cut (set the
  // parameter to NULL if the cost is not desired).
  bool ComputeCut(const std::unordered_map<std::pair<T, T>, double>& edges,
                  std::unordered_set<T>* subgraph1,
                  std::unordered_set<T>* subgraph2,
                  double* cost_or_null) {
    // Create a mapping of node id to index within our linear system.
    IndexNodeIds(edges);

    edge_weight_.resize(node_to_index_map_.size(), node_to_index_map_.size());
    node_weight_.resize(node_to_index_map_.size(), node_to_index_map_.size());
    node_weight_inv_sqrt_.resize(node_to_index_map_.size(),
                                 node_to_index_map_.size());

    // Remove the non-zero entries from any previous calls to this function.
    edge_weight_.setZero();
    node_weight_.setZero();
    node_weight_inv_sqrt_.setZero();

    // Create symmetric weight matrix W where w(i, j) is the weight of the
    // edge
    // between nodes i and j.
    CreateEdgeWeightMatrix(edges);

    // Create diagonal matrix D where d(i) = sum_j w(i, j). Put otherwise, d(i)
    // is the sum of the edge weights connected to node i.
    CreateNodeWeightMatrix();

    // Minimizing the normalized cut is equivalent to finding the vector y
    // such that:
    //
    //   y^t * (D - W) * y
    //   _________________
    //      y^t * D * y
    //
    // is minimized. This is equivalent to a Rayleigh quotient which can be
    // minimized with a generalized eigenvalue system:
    //
    //   (D - W) * y = \lambda * D * y
    //
    // This can be easily transformed into a standard eigenvalue problem:
    //
    //   D^{-1/2} * (D - W) * D^{-1/2} * z = \lambda * z
    //
    // where z = D^{1/2} * y.
    const Eigen::SparseMatrix<double> lhs = node_weight_inv_sqrt_ *
                                            (node_weight_ - edge_weight_) *
                                            node_weight_inv_sqrt_;

    // Note that D^{-1/2} * (D - W) * D^{-1/2} is a symmetric positive
    // semi-definite matrix, so we may use the symmetric eigensolver to find the
    // eigenvalues of lhs.
    SparseSymShiftSolveLLT op(lhs);
    Spectra::SymEigsShiftSolver<double, Spectra::LARGEST_MAGN,
                                SparseSymShiftSolveLLT> eigs(&op, 2, 6, 0.0);
    eigs.init();
    eigs.compute();

    // The eigenvalues will appear in decreasing order. We only care about the
    // eigenvector corresponding to the 2nd smallest eigenvalue.
    const Eigen::VectorXd& z = eigs.eigenvectors().col(0);
    const Eigen::VectorXd y = node_weight_inv_sqrt_ * z;
    FindOptimalCut(y, subgraph1, subgraph2, cost_or_null);

    return true;
  }

 private:
  // Find the optimal cut by searching at evenly spaced points and seeing which
  // point has the lowest cut value. Ideally, the eigenvector y is perfectly
  // split such that the value 0 perfectly divides the graph into the two
  // subgraphs. However, since y was relaxed to be continuous instead of
  // discrete, we need to search for the threshold that splits the eigenvector
  // into the two appropriate groups. We do this by a testing a series of
  // thresholds on the y-vector to determine the grouping and choosing the
  // threshold that produces the lowest cost.
  void FindOptimalCut(const Eigen::VectorXd& y,
                      std::unordered_set<T>* subgraph1,
                      std::unordered_set<T>* subgraph2, double* cost_or_null) {
    const double start_y = y.minCoeff();
    const double stop_y = y.maxCoeff();
    const int num_steps =
        std::min(static_cast<int>(y.size()), options_.num_cuts_to_test);
    const double step_size = (stop_y - start_y) / num_steps;

    double best_cut_value = 0;
    double best_cut_cost = std::numeric_limits<double>::max();

    const auto& node_weight_diag = node_weight_.diagonal();
    const double node_weight_sum = node_weight_diag.sum();
    for (int i = 1; i < num_steps; i++) {
      const double cut_value = start_y + i * step_size;
      // Cut the group based on the cut value such that 1 is in group A and 0 is
      // group B.
      const Eigen::VectorXd cut_grouping =
          (y.array() > cut_value).select(Eigen::VectorXd::Ones(y.size()), 0);

      // Based on our current threshold used for the cut, discretize y so that
      // all values are {1, -b} where
      //   b = \sum_{x_i > 0) d_i / (\sum_{x_i < 0} d_i).
      const double k = node_weight_diag.dot(cut_grouping) / node_weight_sum;
      const double b = k / (1.0 - k);
      const Eigen::VectorXd y_discrete =
          (y.array() > cut_value).select(Eigen::VectorXd::Ones(y.size()), -b);

      // The cost may be computed from y:
      //   ncut cost = y^t * (D - W) * y / (y^t * D * y)
      double cut_cost =
          y_discrete.transpose() * (node_weight_ - edge_weight_) * y_discrete;
      cut_cost /= y_discrete.transpose() * node_weight_ * y_discrete;

      // Select this as the cut if it produces a lower cut cost.
      if (cut_cost < best_cut_cost) {
        best_cut_cost = cut_cost;
        best_cut_value = cut_value;
      }
    }

    // Based on the chosen threshold for the y-values, form the two subgraphs.
    for (const auto& node_id : node_to_index_map_) {
      if (y(node_id.second) <= best_cut_value) {
        subgraph1->emplace(node_id.first);
      } else {
        subgraph2->emplace(node_id.first);
      }
    }

    // Output the cost if desired.
    if (cost_or_null != nullptr) {
      *cost_or_null = best_cut_cost;
    }
  }

  // Create a mapping of node ids to indices that are used for the matrices
  // i.e., which row a particular node id corresponds to.
  void IndexNodeIds(const std::unordered_map<std::pair<T, T>, double>& edges) {
    for (const auto& edge : edges) {
      InsertIfNotPresent(&node_to_index_map_,
                         edge.first.first,
                         node_to_index_map_.size());
      InsertIfNotPresent(&node_to_index_map_,
                         edge.first.second,
                         node_to_index_map_.size());
    }
  }

  // Creates the symmetric edge weight matrix such that
  // w(i,j) = edge_weight(i, j). Only the upper portion of the symmetric
  // matrix is stored.
  void CreateEdgeWeightMatrix(
      const std::unordered_map<std::pair<T, T>, double>& edges) {
    std::vector<Eigen::Triplet<double> > edge_weight_coefficients;
    edge_weight_coefficients.reserve(edges.size());
    for (const auto& edge : edges) {
      int row = FindOrDie(node_to_index_map_, edge.first.first);
      int col = FindOrDie(node_to_index_map_, edge.first.second);
      // Only store the upper part of the matrix.
      if (row > col) {
        std::swap(row, col);
      }

      // Add an entry for w(i, j).
      edge_weight_coefficients.emplace_back(row, col, edge.second);
      edge_weight_coefficients.emplace_back(col, row, edge.second);
    }

    edge_weight_.setFromTriplets(edge_weight_coefficients.begin(),
                                 edge_weight_coefficients.end());
  }

  // Creates the diagonal node weight matrix such that d(i) = sum_j w(i, j)
  void CreateNodeWeightMatrix() {
    std::vector<Eigen::Triplet<double> > node_weight_coefficients;
    std::vector<Eigen::Triplet<double> > node_weight_inv_sqrt_coefficients;
    node_weight_coefficients.reserve(node_to_index_map_.size());
    node_weight_inv_sqrt_coefficients.reserve(node_to_index_map_.size());
    for (const auto& node : node_to_index_map_) {
      // The sum of all edge weights connected to node i is equal to the sum of
      // col(i) in the edge weight matrix.
      const double d_i = edge_weight_.col(node.second).sum();
      node_weight_coefficients.emplace_back(node.first, node.first, d_i);

      // Create D^{-1/2}. Since D is a diagonal matrix the inverse is simply the
      // reciprical of each diagonal matrix.
      node_weight_inv_sqrt_coefficients.emplace_back(node.first, node.first,
                                                     std::sqrt(1.0 / d_i));
    }

    node_weight_.setFromTriplets(node_weight_coefficients.begin(),
                                 node_weight_coefficients.end());
    node_weight_inv_sqrt_.setFromTriplets(
        node_weight_inv_sqrt_coefficients.begin(),
        node_weight_inv_sqrt_coefficients.end());
  }

 private:
  Options options_;
  std::unordered_map<T, int> node_to_index_map_;
  Eigen::SparseMatrix<double> edge_weight_, node_weight_, node_weight_inv_sqrt_;
};

}  // namespace theia

#endif  // THEIA_MATH_GRAPH_NORMALIZED_GRAPH_CUT_H_
