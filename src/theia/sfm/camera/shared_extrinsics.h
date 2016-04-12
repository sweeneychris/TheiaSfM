#ifndef THEIA_CAMERA_SHARED_EXTRINSICS_H
#define THEIA_CAMERA_SHARED_EXTRINSICS_H

#include <Eigen/Core>
#include <Eigen/Geometry>

namespace theia {

class SharedExtrinsics {
public:
  SharedExtrinsics();

  const double* extrinsics() const { return parameters; }
  double* mutable_extrinsics() { return parameters;  }
  enum ExternalParametersIndex {
    POSITION = 0,
    ORIENTATION = 3
  };
  static const int kExtrinsicsSize = 6;
private:
  double parameters[kExtrinsicsSize];
};

}

#endif //THEIA_CAMERA_SHARED_EXTRINSICS_H
