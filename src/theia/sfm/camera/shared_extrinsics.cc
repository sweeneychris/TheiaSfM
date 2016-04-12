#include "shared_extrinsics.h"

namespace theia {

using Eigen::Map;
using Eigen::Matrix;

SharedExtrinsics::SharedExtrinsics() {
  // Set rotation and position to zero (i.e. identity).
  Map<Matrix<double, 1, 6> >(mutable_extrinsics()).setZero();
}

}