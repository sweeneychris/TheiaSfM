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
#include <Eigen/Eigenvalues>
#include <Eigen/SparseCore>
#include <glog/logging.h>

#include <algorithm>
#include <limits>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "cppoptlib/solver/neldermeadsolver.h"
#include "spectra/include/Util/CompInfo.h"
#include "spectra/include/MatOp/SparseCholesky.h"
#include "spectra/include/SymGEigsSolver.h"

#include "theia/math/matrix/spectra_linear_operator.h"
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

    // Remove the non-zero entries from any previous calls to this function.
    edge_weight_.setZero();
    node_weight_.setZero();

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
    // This can be solved using Rayleigh iterations with a sparse eigensolver.
    const Eigen::SparseMatrix<double> lhs = node_weight_ - edge_weight_;
    Spectra::SparseSymMatProd<double> lhs_op(lhs);
    Spectra::SparseCholesky<double> rhs_op(node_weight_);
    Spectra::SymGEigsSolver<double,
                            Spectra::SMALLEST_MAGN,
                            Spectra::SparseSymMatProd<double>,
                            Spectra::SparseCholesky<double>,
                            Spectra::GEIGS_CHOLESKY>
        eigs(&lhs_op, &rhs_op, 2, 6);
    eigs.init();
    eigs.compute();

    // The eigenvalues will appear in decreasing order. We only care about the
    // eigenvector corresponding to the 2nd smallest eigenvalue.
    const Eigen::VectorXd& y = eigs.eigenvectors().col(0);
    FindOptimalCut(y, subgraph1, subgraph2, cost_or_null);

    return true;
  }

 private:
  // The optimal partition point is determined by using a Simplex algorithm to
  // search for the partition that gives the optimal NCut value. This helper
  // class simply computes the NCut cost for the optimization algorithm.
  class ComputeNCutCost : public cppoptlib::Problem<double, 1> {
   public:
    using cppoptlib::Problem<double, 1>::TVector;

    ComputeNCutCost(const Eigen::SparseMatrix<double>& node_weight,
                    const Eigen::SparseMatrix<double>& edge_weight,
                    const Eigen::VectorXd& y)
        : node_weight_(node_weight),
          edge_weight_(edge_weight),
          y_(y),
          node_weight_diag_(node_weight_.diagonal()),
          node_weight_sum_(node_weight_diag_.sum()) {}

    double value(const TVector& cut_value) {
      // Cut the group based on the cut value such that 1 is in group A and 0 is
      // group B.
      const Eigen::VectorXd cut_grouping =
          (y_.array() > cut_value(0))
              .select(Eigen::VectorXd::Ones(y_.size()), 0);

      // Based on our current threshold used for the cut, discretize y so that
      // all values are {1, -b} where
      //   b = \sum_{x_i > 0) d_i / (\sum_{x_i < 0} d_i).
      const double k = node_weight_diag_.dot(cut_grouping) / node_weight_sum_;
      const double b = k / (1.0 - k);
      const Eigen::VectorXd y_discrete =
          (y_.array() > cut_value(0))
              .select(Eigen::VectorXd::Ones(y_.size()), -b);

      // The cost may be computed from y:
      //   ncut cost = y^t * (D - W) * y / (y^t * D * y)
      double cut_cost =
          y_discrete.transpose() * (node_weight_ - edge_weight_) * y_discrete;
      cut_cost /= y_discrete.transpose() * node_weight_ * y_discrete;
      return cut_cost;
    }
   private:
    const Eigen::SparseMatrix<double>& node_weight_;
    const Eigen::SparseMatrix<double>& edge_weight_;
    const Eigen::VectorXd& y_;
    const Eigen::VectorXd node_weight_diag_;
    const double node_weight_sum_;
  };

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
                      std::unordered_set<T>* subgraph2,
                      double* cost_or_null) {
    // We solve for the optimal partition value by performing a minimization
    // using the Nelder-Mead simplex algorithm. The optimal partition will
    // minimize the NCut value (only a local minima is guaranteed).
    ComputeNCutCost ncut_cost(node_weight_, edge_weight_, y);
    // Set the initial value to be the mean of the eigenvector entries.
    typename ComputeNCutCost::TVector best_cut_value(1);
    best_cut_value(0) = y.mean();
    // Try to find partition value where the NCut value is minimized.
    cppoptlib::NelderMeadSolver<ComputeNCutCost> solver;
    solver.minimize(ncut_cost, best_cut_value);

    // Based on the chosen threshold for the y-values, form the two subgraphs.
    for (const auto& node_id : node_to_index_map_) {
      if (y(node_id.second) <= best_cut_value(0)) {
        subgraph1->emplace(node_id.first);
      } else {
        subgraph2->emplace(node_id.first);
      }
    }

    // Output the cost if desired.
    if (cost_or_null != nullptr) {
      *cost_or_null = ncut_cost(best_cut_value);
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
    node_weight_coefficients.reserve(node_to_index_map_.size());
    for (const auto& node : node_to_index_map_) {
      // The sum of all edge weights connected to node i is equal to the sum of
      // col(i) in the edge weight matrix.
      const double d_i = edge_weight_.row(node.second).sum();
      node_weight_coefficients.emplace_back(node.second, node.second, d_i);
    }

    node_weight_.setFromTriplets(node_weight_coefficients.begin(),
                                 node_weight_coefficients.end());
  }

 private:
  Options options_;
  std::unordered_map<T, int> node_to_index_map_;
  Eigen::SparseMatrix<double> edge_weight_, node_weight_;
};

}  // namespace theia

#endif  // THEIA_MATH_GRAPH_NORMALIZED_GRAPH_CUT_H_
