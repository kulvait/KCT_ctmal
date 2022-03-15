#pragma once
// Logging
#include <plog/Log.h>

// External libs
#include <mkl.h>

// Internal libs
#include "GEOMETRY/Geometry3DParallelCameraMatrix.hpp"
#include "MATRIX/utils.hpp"

// Interfaces
#include "GEOMETRY/Geometry3DParallelI.hpp"

namespace KCT {
namespace geometry {
/**
 *Class to represent parallel ray geometry by means of 2x4 homogeneous matrix
 *and detector tilt.
 */
class Geometry3DParallel : public Geometry3DParallelI {
  public:
  Geometry3DParallel(const Geometry3DParallelCameraMatrix &pm,
                     const double cosDetectorTilt);
  // Analogous to ASTRA parallel3d_vec
  Geometry3DParallel(const std::array<double, 3> rayDirection,
                     const std::array<double, 3> detectorOrigin,
                     const std::array<double, 3> VX,
                     const std::array<double, 3> VY);
  // Analogous to ASTRA parallel3d
  Geometry3DParallel(const uint32_t anglesCount, const uint32_t angleNum,
                     const double x_spacing, const double y_spacing);
  bool operator==(const Geometry3DParallel &rhs) {
    if (cosDetectorTilt == rhs.cosDetectorTilt &&
        parallelCameraMatrix == rhs.parallelCameraMatrix) {
      return true;
    } else {
      return false;
    }
  }
  // Implements Geometry3DParallelI
  double pixelSkew() const;
  double detectorTilt() const;
  void directionVectorVR(double *vector3) const;
  void directionVectorVX(double *vector3) const;
  void directionVectorVY(double *vector3) const;
  std::array<double, 3> directionVectorVR() const;
  std::array<double, 3> directionVectorVX() const;
  std::array<double, 3> directionVectorVY() const;
  void normalToDetector(double *vector3) const;
  std::array<double, 3> normalToDetector() const;
  void backprojectToPosition(const double pi, const double pj,
                             double *vector3) const;
  std::array<double, 3> backprojectToPosition(const double pi,
                                              const double pj) const;

  void project(const double x0, const double x1, const double x2, double *pi,
               double *pj) const;
  void project(const float x0, const float x1, const float x2, float *pi,
               float *pj) const;

  void projectionMatrixAsVector8(double *vector8) const;
  std::array<double, 8> projectionMatrixAsVector8() const;

  private:
  Geometry3DParallelCameraMatrix parallelCameraMatrix; // Projection matrix
  double cosDetectorTilt;
  const double zeroPrecisionTolerance = 1e-10;
};

} // namespace matrix
} // namespace KCT
