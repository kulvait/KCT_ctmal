#include "GEOMETRY/Geometry3DParallel.hpp"

namespace KCT {
namespace geometry {
/**
 *Class to represent parallel ray geometry by means of 2x4 homogeneous matrix
 *and detector tilt.
 */
Geometry3DParallel::Geometry3DParallel(
    const Geometry3DParallelCameraMatrix &pcm, const double cosDetectorTilt)
    : parallelCameraMatrix(pcm), cosDetectorTilt(cosDetectorTilt) {}
// Analogous to ASTRA parallel3d_vec
Geometry3DParallel::Geometry3DParallel(
    const std::array<double, 3> rayDirection,
    const std::array<double, 3> detectorOrigin, const std::array<double, 3> VX,
    const std::array<double, 3> VY)
    : parallelCameraMatrix(rayDirection, detectorOrigin, VX, VY) {
  // Compute cosDetectorTilt
  std::array<double, 3> normalToDetector = vectorProduct(VX, VY);
  normalToDetector = normalizeVector<3>(normalToDetector);
  std::array<double, 3> unitRayDirection =
      parallelCameraMatrix.directionVectorVR();
  cosDetectorTilt =
      std::abs(vectorDotProduct<3>(normalToDetector, unitRayDirection));
}
// Analogous to ASTRA parallel3d
Geometry3DParallel::Geometry3DParallel(const uint32_t anglesCount,
                                       const uint32_t angleNum,
                                       const double x_spacing,
                                       const double y_spacing)
    : parallelCameraMatrix(anglesCount, angleNum, x_spacing, y_spacing) {
  cosDetectorTilt = 1.0;
}

// Implements Geometry3DParallelI
double Geometry3DParallel::pixelSkew() const {
  // I have to normalize VX and VY and compute their scalar product
  std::array<double, 3> VXN =
      normalizeVector<3>(parallelCameraMatrix.directionVectorVX());
  std::array<double, 3> VYN =
      normalizeVector<3>(parallelCameraMatrix.directionVectorVY());
  // Here without absolute value as the vectors VX and VY direction matters
  return vectorDotProduct<3>(VXN, VYN);
}
double Geometry3DParallel::detectorTilt() const { return cosDetectorTilt; }
void Geometry3DParallel::directionVectorVR(double *vector3) const {
  parallelCameraMatrix.directionVectorVR(vector3);
}
void Geometry3DParallel::directionVectorVX(double *vector3) const {
  parallelCameraMatrix.directionVectorVX(vector3);
}
void Geometry3DParallel::directionVectorVY(double *vector3) const {
  parallelCameraMatrix.directionVectorVY(vector3);
}
/**
*
* @return Unit vector in the direction of incomming rays or its oposite.
*/
std::array<double, 3> Geometry3DParallel::directionVectorVR() const {
  return parallelCameraMatrix.directionVectorVR();
}
std::array<double, 3> Geometry3DParallel::directionVectorVX() const {
  return parallelCameraMatrix.directionVectorVX();
}
std::array<double, 3> Geometry3DParallel::directionVectorVY() const {
  return parallelCameraMatrix.directionVectorVY();
}
void Geometry3DParallel::backprojectToPosition(const double PX, const double PY,
                                               double *vector3) const {
  parallelCameraMatrix.backprojectToPosition(PX, PY, vector3);
}

std::array<double, 3>
Geometry3DParallel::backprojectToPosition(const double PX,
                                          const double PY) const {
  return parallelCameraMatrix.backprojectToPosition(PX, PY);
}

void Geometry3DParallel::project(const double x0, const double x1,
                                 const double x2, double *PX,
                                 double *PY) const {
  parallelCameraMatrix.project(x0, x1, x2, PX, PY);
}
void Geometry3DParallel::project(const float x0, const float x1, const float x2,
                                 float *PX, float *PY) const {
  parallelCameraMatrix.project(x0, x1, x2, PX, PY);
}

void Geometry3DParallel::projectionMatrixAsVector8(double *vector8) const {
  parallelCameraMatrix.projectionMatrixAsVector8(vector8);
}

std::array<double, 8> Geometry3DParallel::projectionMatrixAsVector8() const {
  return parallelCameraMatrix.projectionMatrixAsVector8();
}

} // namespace matrix
} // namespace KCT
