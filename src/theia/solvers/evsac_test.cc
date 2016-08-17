// Copyright (C) 2014 The Regents of the University of California (Regents).
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
// Author: Victor Fragoso (vfragoso@cs.ucsb.edu)

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <algorithm>
#include <utility>
#include <vector>

#include "gtest/gtest.h"
#include "theia/solvers/estimator.h"
#include "theia/solvers/evsac.h"
#include "theia/sfm/pose/four_point_homography.h"
#include "theia/util/random.h"
#include "theia/test/test_utils.h"

namespace theia {
using std::vector;
using std::pair;
using Eigen::Vector2d;
using Eigen::RowVectorXd;
using Eigen::RowVector4d;
using Eigen::Vector3d;
using Eigen::MatrixXd;
using Eigen::Matrix3d;

namespace {

RandomNumberGenerator rng(60);

// Rayleigh samples taken from MATLAB with sigma=1.
const vector<double> rayl_samples = {
  9.442826e-01, 1.571453e+00, 9.874990e-01, 1.803729e+00, 7.826802e-01,
  8.438893e-01, 1.877576e-01, 2.045208e+00, 9.909384e-01, 8.541785e-01,
  1.048271e+00, 1.136111e+00, 7.069109e-01, 3.450315e-01, 1.019216e+00,
  1.069638e+00, 8.893999e-01, 1.829623e+00, 1.235894e+00, 2.195386e-01,
  2.303448e+00, 9.455074e-01, 9.461517e-01, 1.795313e+00, 1.310188e+00,
  9.176731e-01, 3.815293e-01, 1.213376e+00, 8.197362e-01, 3.146557e-01
};

// mixture samples: matrix of 50 x 6 containing in the first column samples of a
// gamma(3, 2) and the 2:6 columns are sorted samples taken from a chi2rnd(128).
// This is a separable case. The data is written using a column major.
// This data was generated with MATLAB.
const vector<double> mixture_samples = {
  8.453391e+00, 6.643481e+00, 7.073280e+00, 7.620230e-01, 4.796192e+00,
  4.265579e+00, 7.504951e+00, 6.491244e+00, 1.086329e+01, 1.008544e+01,
  4.580231e+00, 2.009185e+00, 4.833628e+00, 2.524501e+00, 5.473356e+00,
  1.412114e+01, 7.078814e+00, 5.815086e+00, 5.774830e+00, 7.891432e+00,
  4.030500e+00, 3.868219e+00, 7.489940e+00, 1.655775e+00, 2.400819e+00,
  3.553688e+00, 1.342447e+01, 5.690786e+00, 1.377225e+01, 2.944703e+00,
  6.516840e+00, 4.445703e+00, 8.211898e+00, 3.078726e+00, 1.199771e+01,
  4.297356e+00, 3.419880e+00, 3.356973e+00, 1.402518e+01, 6.971956e+00,
  4.778875e+00, 5.016027e+00, 1.416609e+00, 4.515507e+00, 6.002494e+00,
  1.032861e+01, 5.377348e+00, 5.531420e+00, 5.714113e+00, 2.514066e+00,
  7.747292e+01, 8.828736e+01, 7.707717e+01, 6.252057e+01, 8.828689e+01,
  7.832088e+01, 8.107225e+01, 8.690371e+01, 8.368767e+01, 8.400832e+01,
  7.684186e+01, 7.281568e+01, 8.273788e+01, 8.661383e+01, 7.883183e+01,
  8.208945e+01, 7.917460e+01, 7.973665e+01, 8.819710e+01, 8.426296e+01,
  8.509478e+01, 8.832398e+01, 8.264218e+01, 8.444737e+01, 7.573799e+01,
  8.841938e+01, 7.811256e+01, 7.576132e+01, 8.681609e+01, 8.489867e+01,
  8.107643e+01, 8.662812e+01, 8.254524e+01, 7.886447e+01, 8.168945e+01,
  7.817184e+01, 8.283279e+01, 8.639165e+01, 7.494652e+01, 9.019843e+01,
  8.972476e+01, 8.043246e+01, 7.397081e+01, 8.708566e+01, 8.882751e+01,
  8.366903e+01, 7.833746e+01, 8.907670e+01, 8.344260e+01, 8.655409e+01,
  8.729785e+01, 8.922214e+01, 7.845613e+01, 7.889089e+01, 9.060706e+01,
  8.381995e+01, 8.526234e+01, 8.729997e+01, 8.701973e+01, 8.490448e+01,
  8.105177e+01, 8.261784e+01, 8.548979e+01, 8.890356e+01, 8.402725e+01,
  8.315342e+01, 8.384340e+01, 8.553590e+01, 8.991483e+01, 8.788811e+01,
  8.534967e+01, 8.863323e+01, 8.544508e+01, 8.514051e+01, 8.376380e+01,
  8.991751e+01, 8.346599e+01, 8.031017e+01, 9.317175e+01, 8.753779e+01,
  8.496033e+01, 8.801973e+01, 8.618550e+01, 7.970668e+01, 8.680145e+01,
  8.072150e+01, 8.624960e+01, 8.651754e+01, 8.669568e+01, 9.083972e+01,
  8.974145e+01, 8.887377e+01, 8.312898e+01, 9.055060e+01, 8.900577e+01,
  8.639638e+01, 7.949458e+01, 9.038057e+01, 8.498672e+01, 8.846037e+01,
  8.816720e+01, 9.004428e+01, 8.555926e+01, 8.459183e+01, 9.204340e+01,
  8.459763e+01, 8.829754e+01, 8.779764e+01, 8.741255e+01, 8.579928e+01,
  8.524173e+01, 9.019785e+01, 8.968418e+01, 9.033841e+01, 8.494151e+01,
  8.476563e+01, 8.703247e+01, 8.642406e+01, 9.083866e+01, 8.802507e+01,
  8.853094e+01, 9.004607e+01, 8.565775e+01, 8.636985e+01, 8.489890e+01,
  9.017141e+01, 8.671391e+01, 8.102536e+01, 9.395337e+01, 8.863653e+01,
  8.718898e+01, 9.115773e+01, 8.737245e+01, 8.339203e+01, 8.775825e+01,
  8.745352e+01, 8.636032e+01, 8.897125e+01, 8.724385e+01, 9.114361e+01,
  9.073200e+01, 8.902639e+01, 8.513250e+01, 9.055553e+01, 9.092100e+01,
  8.699845e+01, 8.183630e+01, 9.096764e+01, 8.619657e+01, 8.918183e+01,
  8.895334e+01, 9.072155e+01, 8.596572e+01, 8.686888e+01, 9.209827e+01,
  8.697255e+01, 8.974051e+01, 8.995643e+01, 9.067128e+01, 8.862904e+01,
  8.845180e+01, 9.219812e+01, 9.056826e+01, 9.044349e+01, 8.776265e+01,
  8.749650e+01, 8.801521e+01, 8.643755e+01, 9.125138e+01, 9.063493e+01,
  8.882734e+01, 9.186411e+01, 8.697559e+01, 8.640061e+01, 8.496423e+01,
  9.220204e+01, 8.738610e+01, 8.893682e+01, 9.421022e+01, 9.167722e+01,
  8.816749e+01, 9.155549e+01, 8.850337e+01, 8.341849e+01, 8.962811e+01,
  9.011022e+01, 8.938827e+01, 8.948222e+01, 8.724831e+01, 9.150108e+01,
  9.092604e+01, 8.930900e+01, 8.946705e+01, 9.060121e+01, 9.120959e+01,
  8.773249e+01, 8.563068e+01, 9.159597e+01, 8.808799e+01, 8.918286e+01,
  8.926355e+01, 9.278073e+01, 9.017833e+01, 8.710429e+01, 9.319228e+01,
  8.887087e+01, 8.999791e+01, 9.042229e+01, 9.158277e+01, 9.135374e+01,
  8.878243e+01, 9.259866e+01, 9.157813e+01, 9.142548e+01, 8.943605e+01,
  8.861633e+01, 8.915538e+01, 8.646916e+01, 9.144094e+01, 9.096172e+01,
  8.895131e+01, 9.328853e+01, 8.875806e+01, 8.644363e+01, 8.676286e+01,
  9.310095e+01, 8.763685e+01, 8.922829e+01, 9.451073e+01, 9.261248e+01,
  8.873643e+01, 9.319027e+01, 8.863945e+01, 8.480577e+01, 9.076276e+01,
  9.062156e+01, 8.950944e+01, 9.012737e+01, 9.051290e+01, 9.264635e+01,
  9.108225e+01, 9.045001e+01, 9.012979e+01, 9.157871e+01, 9.429332e+01,
  8.795334e+01, 8.852589e+01, 9.203239e+01, 8.977151e+01, 8.969222e+01
};

// Homography to be estimated; in a column major order.
const vector<double> homography_values = {
  6.682947e-01, -4.902747e-02, -1.899965e-04, 1.834432e-02, 8.935012e-01,
  2.069200e-06, 3.951916e+01, 6.181007e+01, 1
};

// Correspondence matrix 100 x 4 in a column major. The first block of 50x4 are
// correct correspondences, the rest are outliers. These coordinates have
// independent Gaussian noise with a std. deviation of 1.0.
const vector<double> correspondences_vec = {
  2.109930e+02, 4.929580e+02, 9.219877e+00, 2.498119e+01, 1.002647e+02,
  3.898529e+02, 4.390821e+02, 3.870780e+02, 2.704160e+02, 3.266696e+02,
  1.766281e+02, 4.478900e+02, 1.125440e+02, 4.125498e+02, 1.105928e+02,
  2.226137e+02, 3.747635e+02, 4.680257e+02, 4.851485e+01, 5.566712e+02,
  4.665804e+02, 2.893302e+02, 2.600961e+02, 2.688608e+02, 1.836158e+02,
  3.032625e+02, 3.059163e+02, 4.923486e+02, 4.775595e+02, 3.876093e+02,
  2.271383e+02, 4.851530e+02, 3.199514e+02, 2.115324e+02, 5.634389e+02,
  5.247111e+02, 3.286959e+02, 3.724085e+02, 3.509394e+02, 1.252495e+02,
  1.813015e+02, 2.824793e+02, 1.384568e+02, 5.062508e+02, 1.159545e+02,
  1.359150e+02, 1.020835e+02, 1.382683e+02, 2.604055e+02, 1.861664e+02,
  4.645975e+00, 2.903209e+01, 4.000611e+02, 3.617305e+02, 3.149728e+02,
  4.358999e+02, 4.239178e+02, 4.683563e+02, 1.721794e+02, 4.149753e+02,
  3.319141e+02, 2.393220e+02, 3.669097e+01, 4.670799e+02, 2.016830e+02,
  3.667210e+02, 4.440568e+02, 6.210350e+01, 7.524254e+01, 3.296017e+02,
  2.889616e+02, 5.344890e+02, 4.804038e+02, 4.409107e+02, 3.013919e+01,
  4.272096e+01, 5.539807e+01, 4.797210e+02, 5.648239e+02, 4.109137e+02,
  7.850858e+01, 4.345643e+02, 6.581732e+01, 7.139581e+01, 3.836461e+02,
  1.966248e+02, 3.930609e+02, 4.483654e+02, 3.497904e+02, 4.423594e+02,
  1.422316e+02, 4.397186e+02, 5.800697e+02, 5.190700e+02, 5.203616e+01,
  2.211205e+02, 2.218160e+02, 4.113084e+02, 3.601563e+02, 4.748742e+02,
  5.542944e+02, 2.588970e+02, 1.120773e+02, 5.428586e+02, 5.878013e+02,
  2.630485e+02, 6.619341e+01, 1.556820e+02, 2.444180e+02, 3.572979e+02,
  1.574801e+02, 3.610997e+02, 4.259294e+02, 1.309186e+02, 6.983930e+01,
  1.770379e+02, 1.911712e+02, 2.544043e+02, 3.031367e+02, 5.070228e+01,
  1.584842e+02, 4.796945e+02, 1.964769e+01, 5.564088e+02, 4.385060e+02,
  2.937751e+02, 3.450911e+02, 1.405259e+02, 2.760759e+02, 5.766825e+02,
  3.274650e+02, 3.125545e+02, 1.387690e+02, 2.941405e+02, 3.742952e+02,
  4.065894e+02, 2.379792e+02, 2.198932e+02, 5.926096e+02, 2.353479e+01,
  5.298896e+02, 5.490397e+02, 4.788007e+02, 5.832287e+01, 1.568323e+02,
  2.013969e+02, 4.065614e+02, 8.251036e+01, 4.338126e+02, 6.559317e+01,
  2.196165e+02, 1.253937e+02, 5.081987e+01, 4.633945e+02, 1.238142e+02,
  2.326183e+02, 3.291930e+02, 1.381996e+02, 3.845032e+02, 2.894556e+02,
  9.051599e+01, 4.683226e+02, 5.972502e+01, 1.758393e+02, 1.422009e+02,
  3.202918e+02, 5.458266e+01, 2.416103e+02, 6.482667e+01, 6.758079e+01,
  4.712534e+02, 1.743481e+02, 3.615681e+02, 5.796173e+02, 2.581659e+02,
  4.155404e+02, 4.538116e+02, 2.597600e+02, 3.932358e+02, 6.598496e+01,
  5.612804e+02, 1.124848e+02, 1.601972e+02, 4.778079e+02, 2.929455e+02,
  4.618749e+02, 2.381374e+02, 1.641764e+02, 2.213293e+01, 4.024450e+02,
  2.583774e+02, 2.713724e+02, 3.670860e+02, 3.595057e+01, 1.914731e+02,
  4.627083e+02, 4.179141e+02, 7.392978e+01, 7.760363e+01, 5.605467e+01,
  1.987635e+02, 4.112680e+02, 4.912180e+01, 6.735509e+01, 1.199456e+02,
  3.309524e+02, 3.641270e+02, 3.265804e+02, 2.366130e+02, 2.831971e+02,
  1.675177e+02, 3.760588e+02, 1.250890e+02, 3.456811e+02, 1.168375e+02,
  1.978597e+02, 3.163277e+02, 3.906989e+02, 7.774040e+01, 4.616683e+02,
  3.890552e+02, 2.577773e+02, 2.261382e+02, 2.420492e+02, 1.761295e+02,
  2.619288e+02, 2.651394e+02, 4.092051e+02, 3.999343e+02, 3.323751e+02,
  2.066042e+02, 4.077947e+02, 2.710809e+02, 1.933171e+02, 4.738081e+02,
  4.410128e+02, 2.807206e+02, 3.150165e+02, 3.039743e+02, 1.247797e+02,
  1.760791e+02, 2.503036e+02, 1.435406e+02, 4.186917e+02, 1.204414e+02,
  1.360157e+02, 1.173695e+02, 1.348523e+02, 2.328355e+02, 1.717321e+02,
  5.889017e+00, 2.548464e+02, 3.929143e+02, 4.344113e+02, 3.169716e+02,
  6.495806e+01, 3.789487e+02, 7.373572e+01, 7.929252e+01, 5.702445e+01,
  8.626114e+01, 1.021039e+02, 1.193326e+02, 1.901230e+02, 1.899052e+02,
  1.316620e+02, 1.507456e+02, 5.378136e+02, 4.221004e+02, 3.347878e+02,
  1.109960e+02, 1.278485e+02, 4.770814e+01, 5.496866e+02, 4.236471e+02,
  3.357621e+02, 1.880784e+02, 9.946568e+01, 3.716255e+02, 5.935969e+02,
  1.024503e+02, 1.555056e+02, 2.386835e+02, 4.403885e+01, 4.095528e+02,
  2.418289e+02, 5.898439e+02, 2.426534e+02, 3.726701e+02, 9.310415e+01,
  2.279604e+02, 9.561463e+01, 4.550931e+02, 5.227041e+02, 2.082805e+02,
  4.120454e+02, 1.754124e+02, 3.183182e+02, 4.978882e+02, 3.584427e+02,
  5.678280e+02, 2.957596e+02, 1.607652e+02, 5.465260e+02, 5.931582e+02,
  2.999853e+02, 1.107222e+02, 1.946621e+02, 2.817496e+02, 3.881044e+02,
  1.999151e+02, 3.977503e+02, 4.478460e+02, 1.740162e+02, 1.204976e+02,
  2.195291e+02, 2.299141e+02, 2.920531e+02, 3.345095e+02, 9.109464e+01,
  1.981237e+02, 5.032883e+02, 6.603967e+01, 5.751236e+02, 4.587299e+02,
  3.280750e+02, 3.794400e+02, 1.812523e+02, 3.137526e+02, 6.025345e+02,
  3.607918e+02, 3.493377e+02, 1.821647e+02, 3.272579e+02, 4.137077e+02,
  4.445427e+02, 2.745654e+02, 2.592997e+02, 6.142718e+02, 7.956021e+01,
  5.451015e+02, 5.663842e+02, 4.946384e+02, 9.817405e+01, 2.019947e+02,
  2.400915e+02, 4.303666e+02, 1.309441e+02, 4.566980e+02, 1.144410e+02,
  2.002009e+02, 1.789359e+02, 2.706043e+02, 2.525987e+02, 2.144653e+02,
  3.344873e+02, 4.467363e+02, 2.551269e+02, 2.583241e+02, 7.490351e+01,
  1.717659e+01, 1.734226e+02, 1.917575e+02, 3.925656e+02, 5.739818e+02,
  5.627405e+02, 2.736563e+02, 1.437468e+02, 4.556250e+02, 4.557298e+02,
  4.446600e+02, 4.476361e+02, 6.481773e+01, 4.090091e+02, 2.777925e+02,
  1.280586e+02, 5.900945e+01, 4.945155e+02, 1.030594e+02, 9.734951e+01,
  4.001592e+02, 5.369335e+02, 3.110142e+02, 4.202034e+02, 9.131150e+01,
  5.719054e+02, 3.225368e+02, 4.081865e+02, 2.081653e+01, 4.868790e+02,
  4.495029e+02, 7.125297e+01, 3.158927e+02, 1.946031e+02, 3.265794e+02,
  2.391768e+02, 2.513276e+02, 1.089671e+02, 1.523320e+02, 1.104546e+01
};

// Data generated to complement the vector of mixture_samples used for testing
// the estimation of homography.
const vector<double> outlier_distances = {
  8.895148e+01, 7.953063e+01, 8.460098e+01, 8.994936e+01, 8.566362e+01,
  8.242079e+01, 8.163974e+01, 7.999375e+01, 7.824075e+01, 8.478654e+01,
  8.227023e+01, 8.435697e+01, 8.541193e+01, 7.860078e+01, 8.512051e+01,
  7.880110e+01, 8.307795e+01, 8.658602e+01, 7.996277e+01, 8.415694e+01,
  8.123535e+01, 7.934643e+01, 8.594829e+01, 8.829724e+01, 8.380728e+01,
  7.656431e+01, 8.495329e+01, 8.480894e+01, 8.343034e+01, 8.524496e+01,
  8.385521e+01, 8.348227e+01, 8.162923e+01, 8.547007e+01, 8.363076e+01,
  8.392757e+01, 8.468554e+01, 7.392006e+01, 8.706441e+01, 8.307758e+01,
  8.729904e+01, 8.486077e+01, 8.624949e+01, 7.646400e+01, 8.077639e+01,
  8.789022e+01, 7.669106e+01, 8.373965e+01, 8.512865e+01, 9.049012e+01,
  9.086200e+01, 8.598175e+01, 8.654790e+01, 9.051895e+01, 8.810061e+01,
  8.440460e+01, 8.803212e+01, 8.244076e+01, 8.480498e+01, 8.774113e+01,
  8.609864e+01, 9.058883e+01, 8.547178e+01, 8.423424e+01, 8.668341e+01,
  8.304153e+01, 8.406133e+01, 9.046712e+01, 8.601846e+01, 8.436973e+01,
  8.333788e+01, 8.287855e+01, 8.640524e+01, 8.913714e+01, 8.513370e+01,
  8.466741e+01, 8.772383e+01, 8.867448e+01, 8.806910e+01, 8.640077e+01,
  8.732349e+01, 8.568596e+01, 8.481050e+01, 8.773720e+01, 8.398822e+01,
  8.751314e+01, 8.992481e+01, 8.846628e+01, 8.778357e+01, 8.658942e+01,
  9.178961e+01, 8.691123e+01, 9.099425e+01, 8.432981e+01, 8.716866e+01,
  8.862666e+01, 7.767600e+01, 8.544150e+01, 9.084609e+01, 9.129871e+01,
  9.101739e+01, 8.645098e+01, 8.717636e+01, 9.095094e+01, 9.083382e+01,
  8.479854e+01, 9.068514e+01, 9.110348e+01, 8.742729e+01, 8.811386e+01,
  8.997135e+01, 9.062248e+01, 8.964491e+01, 8.957363e+01, 8.845934e+01,
  8.408958e+01, 8.916568e+01, 9.142303e+01, 9.157848e+01, 8.706515e+01,
  8.511824e+01, 8.887819e+01, 8.665248e+01, 9.134160e+01, 8.601966e+01,
  8.592578e+01, 8.949893e+01, 8.915504e+01, 8.982826e+01, 9.063684e+01,
  8.798646e+01, 8.685137e+01, 8.499423e+01, 8.850391e+01, 8.480393e+01,
  8.787008e+01, 9.257978e+01, 9.121866e+01, 8.897564e+01, 8.897865e+01,
  9.310619e+01, 9.036521e+01, 9.102717e+01, 8.526297e+01, 8.844707e+01,
  8.881979e+01, 8.438537e+01, 8.973968e+01, 9.154571e+01, 9.229044e+01,
  9.103770e+01, 8.855978e+01, 8.781302e+01, 9.212046e+01, 9.090739e+01,
  8.499870e+01, 9.181944e+01, 9.274805e+01, 8.749219e+01, 9.011331e+01,
  9.031167e+01, 9.076769e+01, 8.981600e+01, 8.962728e+01, 8.981828e+01,
  8.698774e+01, 8.917790e+01, 9.228209e+01, 9.227437e+01, 8.978103e+01,
  8.786538e+01, 8.916256e+01, 8.836609e+01, 9.144562e+01, 8.611135e+01,
  8.782449e+01, 8.950367e+01, 8.920596e+01, 9.011024e+01, 9.089128e+01,
  8.822924e+01, 8.898492e+01, 8.963555e+01, 9.134662e+01, 8.673048e+01,
  9.189652e+01, 9.359115e+01, 9.142325e+01, 9.175463e+01, 9.053526e+01,
  9.419661e+01, 9.255030e+01, 9.251931e+01, 8.559187e+01, 8.901613e+01,
  8.902638e+01, 8.632767e+01, 9.070651e+01, 9.197341e+01, 9.320143e+01,
  9.204102e+01, 8.897606e+01, 8.962652e+01, 9.278509e+01, 9.127775e+01,
  9.347129e+01, 9.300554e+01, 9.443262e+01, 8.861573e+01, 9.218864e+01,
  9.153521e+01, 9.113444e+01, 8.988825e+01, 8.965697e+01, 9.224423e+01,
  8.841865e+01, 9.239366e+01, 9.235699e+01, 9.281394e+01, 9.113598e+01,
  8.797721e+01, 9.036616e+01, 8.909086e+01, 9.161517e+01, 8.713308e+01,
  8.803306e+01, 9.092268e+01, 9.000626e+01, 9.147762e+01, 9.326346e+01,
  8.979955e+01, 8.941126e+01, 9.192247e+01, 9.216250e+01, 8.856981e+01,
  9.206467e+01, 9.467586e+01, 9.187157e+01, 9.223276e+01, 9.163224e+01,
  9.557799e+01, 9.278091e+01, 9.276446e+01, 8.860777e+01, 8.919436e+01,
  9.001741e+01, 8.785049e+01, 9.162323e+01, 9.229596e+01, 9.358814e+01,
  9.218509e+01, 8.936480e+01, 8.986928e+01, 9.411375e+01, 9.132268e+01,
  9.356743e+01, 9.311709e+01, 9.457162e+01, 9.054907e+01, 9.276146e+01,
  9.181887e+01, 9.141717e+01, 9.091806e+01, 8.969953e+01, 9.288332e+01,
  9.054916e+01, 9.243261e+01, 9.242040e+01, 9.304591e+01, 9.350347e+01,
  8.855729e+01, 9.067934e+01, 9.002463e+01, 9.187813e+01, 8.854240e+01,
  8.944499e+01, 9.218716e+01, 9.051515e+01, 9.208943e+01, 9.344828e+01,
  9.103225e+01, 9.146842e+01, 9.300111e+01, 9.257911e+01, 8.932510e+01,
  9.219331e+01, 9.502834e+01, 9.239294e+01, 9.276069e+01, 9.226127e+01,
  9.566503e+01, 9.321390e+01, 9.283700e+01, 8.981228e+01, 8.927680e+01,
  9.063661e+01, 8.865151e+01, 9.180143e+01, 9.470521e+01, 9.416588e+01
};

typedef pair<Vector2d, Vector2d> Correspondence;

class HomographyEstimator : public Estimator<Correspondence, Matrix3d> {
 public:
  HomographyEstimator() {}
  ~HomographyEstimator() {}

  double SampleSize() const override { return 4; }

  bool EstimateModel(const vector<Correspondence>& data,
                     vector<Matrix3d>* models) const override {
    vector<Vector2d> image_1_points;
    vector<Vector2d> image_2_points;
    for (const Correspondence& correspondence : data) {
      image_1_points.push_back(correspondence.first);
      image_2_points.push_back(correspondence.second);
    }
    Matrix3d homography;
    const bool homography_flag =
        FourPointHomography(image_1_points, image_2_points, &homography);
    homography /= homography(2, 2);
    models->emplace_back(homography);
    return homography_flag;
  }

  double Error(const Correspondence& correspondence,
               const Matrix3d& homography) const override {
    Vector3d x_hat = homography * correspondence.first.homogeneous();
    x_hat /= x_hat(2);
    return (x_hat - correspondence.second.homogeneous()).norm();
  }
};
}  // namespace

// Tests the calculation of the k, the parameter to select k smallest distances.
TEST(EvsacSamplerTest, CalculateKSmallestDistances) {
  // 0.01% of 100 features should return 1!
  const int kNumReferenceFeatures = 100;
  const double kPercentile = 0.01;
  EXPECT_EQ(EvsacSampler<Vector2d>::CalculateKSmallestDistances(
      kNumReferenceFeatures, kPercentile), 1);
}

// Tests the prediction using MRRayleigh.
TEST(EvsacSamplerTest, MRRayleighPrediction) {
  const double kPredictorThreshold = 0.65;
  RowVectorXd correct_vec(rayl_samples.size() + 1);
  RowVectorXd incorrect_vec(rayl_samples.size() + 1);

  // Correct prediction case.
  double* data = correct_vec.data();
  data[0] = 0.01;
  std::copy(rayl_samples.begin(), rayl_samples.end(), data + 1);
  EXPECT_TRUE(
      EvsacSampler<Vector2d>::MRRayleigh(correct_vec, kPredictorThreshold));

  // Incorrect prediction case.
  data = incorrect_vec.data();
  data[0] = 1.00;
  std::copy(rayl_samples.begin(), rayl_samples.end(), data + 1);
  EXPECT_FALSE(
      EvsacSampler<Vector2d>::MRRayleigh(incorrect_vec, kPredictorThreshold));
}

// Tests the calculation of the mixture model parameters.
TEST(EvsacSamplerTest, MixtureModelParamsCalculation) {
  const double kPredictorThreshold = 0.65;
  // Parameters obtained with MATLAB
  // Gamma parameters gamp = [3.2807 1.8495].
  // Reversed GEV parameters gev = [0.0148 4.1872 -84.7617].
  EvsacSampler<Vector2d>::MixtureModelParams mixture_model_params;
  vector<float> probabilities;
  vector<float> sampling_weights;
  const MatrixXd distances =
      Eigen::Map<const MatrixXd>(&mixture_samples[0], 50, 6);
  EvsacSampler<Vector2d>::CalculateMixtureModel(
      distances, kPredictorThreshold, MLE,
      &mixture_model_params, &probabilities, &sampling_weights);
  // Checking mixture model params.
  // Checking Gamma parameters.
  EXPECT_NEAR(mixture_model_params.k, 3.2807, 0.1);
  EXPECT_NEAR(mixture_model_params.theta, 1.8495, 0.1);
  // Checking reversed GEV parameters.
  EXPECT_NEAR(mixture_model_params.xi, 0.0148, 0.1);
  EXPECT_NEAR(mixture_model_params.sigma, 4.1872, 0.1);
  EXPECT_NEAR(mixture_model_params.mu, -84.7617, 0.1);
  // Checking inlier ratio.
  EXPECT_NEAR(mixture_model_params.inlier_ratio, 1.0, 0.1);
  // Check that sampling weights are high and non zero!
  float avg_p = 0.0f;
  for (const float p : sampling_weights) avg_p += p;
  avg_p /= sampling_weights.size();
  EXPECT_GT(avg_p, 0.95);
  VLOG(1) << "Average probability: " << avg_p;
}

// Tests the estimation of Homography using Evsac.
TEST(EvsacTest, HomographyEstimationWithEvsacSampling) {
  vector<double> distance_vector;
  // Copy the distances of the inliers.
  distance_vector.insert(
      distance_vector.end(), mixture_samples.begin(), mixture_samples.end());
  // Copy the distances of the outliers.
  distance_vector.insert(
      distance_vector.end(), outlier_distances.begin(),
      outlier_distances.end());
  const MatrixXd distances =
      Eigen::Map<const MatrixXd>(&distance_vector[0], 100, 6);
  // Load ground truth homography.
  const Matrix3d homography_gt =
      Eigen::Map<const Matrix3d>(&homography_values[0]);
  // Load correspondences.
  const MatrixXd correspondences =
      Eigen::Map<const MatrixXd>(&correspondences_vec[0], 100, 4);
  // Generate the correspondence vector.
  vector<Correspondence> data;
  for (int row = 0; row < correspondences.rows(); ++row) {
    const RowVector4d& row_vec = correspondences.row(row);
    data.emplace_back(
        Vector2d(row_vec(0), row_vec(1)), Vector2d(row_vec(2), row_vec(3)));
  }
  HomographyEstimator homography_estimator;
  Matrix3d homography;
  RansacParameters params;
  params.rng = std::make_shared<RandomNumberGenerator>(rng);
  params.error_thresh = 5.0;  // 5px of error
  const double kPredictorThreshold = 0.65;
  Evsac<HomographyEstimator> evsac_homography(
      params, homography_estimator, distances, kPredictorThreshold, MLE);
  evsac_homography.Initialize();
  RansacSummary summary;
  evsac_homography.Estimate(data, &homography, &summary);
  EXPECT_LT((homography_gt - homography).norm() / 9.0, 1.0);
  VLOG(1) << "Estimate: \n" << homography << "\n GroundTruth: \n"
          << homography_gt;
}
}  // namespace theia
