#include "GEOMETRY/Geometry3DParallelCameraMatrix.hpp"

namespace KCT {
namespace geometry {

void Geometry3DParallelCameraMatrix::computeRayDirection();

Geometry3DParallelCameraMatrix::Geometry3DParallelCameraMatrix(
    const double *vector8);
Geometry3DParallelCameraMatrix::Geometry3DParallelCameraMatrix(
    const double (&initArray)[2 * 4]);
/**Constructor from the 2x4 Matrix.*/
Geometry3DParallelCameraMatrix::Geometry3DParallelCameraMatrix(
    const Matrix &pm);
// ASTRA like initialization
Geometry3DParallelCameraMatrix::Geometry3DParallelCameraMatrix(
    const std::array<double, 3> rayDirection,
    const std::array<double, 3> detectorOrigin, const std::array<double, 3> VX,
    const std::array<double, 3> VY);
Geometry3DParallelCameraMatrix::Geometry3DParallelCameraMatrix(
    const uint32_t anglesCount, const uint32_t angleNum, const double x_spacing,
    const double y_spacing);
bool Geometry3DParallelCameraMatrix::
operator==(const Geometry3DParallelCameraMatrix &rhs) {
  for (int i = 0; i != 2; i++) {
    for (int j = 0; j != 4; j++) {
      if ((*this)(i, j) != (rhs)(i, j))
        return false;
    }
  }
  return true;
}
// Implements Geometry3DParallelCameraMatrixI
double Geometry3DParallelCameraMatrix::pixelSkew() const;
void Geometry3DParallelCameraMatrix::directionVectorVR(double *vector3) const;
void Geometry3DParallelCameraMatrix::directionVectorVX(double *vector3) const;
void Geometry3DParallelCameraMatrix::directionVectorVY(double *vector3) const;
std::array<double, 3> Geometry3DParallelCameraMatrix::directionVectorVR() const;
std::array<double, 3> Geometry3DParallelCameraMatrix::directionVectorVX() const;
std::array<double, 3> Geometry3DParallelCameraMatrix::directionVectorVY() const;
void Geometry3DParallelCameraMatrix::backprojectToPosition(
    const double pi, const double pj, double *vector3) const;
std::array<double, 3>
Geometry3DParallelCameraMatrix::backprojectToPosition(const double pi,
                                                      const double pj) const;
void Geometry3DParallelCameraMatrix::project(const double x0, const double x1,
                                             const double x2, double *PX,
                                             double *PY) const;
void Geometry3DParallelCameraMatrix::project(const float x0, const float x1,
                                             const float x2, float *PX,
                                             float *PY) const;

void Geometry3DParallelCameraMatrix::projectionMatrixAsVector8(double *vector8)
    const;
std::array<double, 8>
Geometry3DParallelCameraMatrix::projectionMatrixAsVector8() const;
std::string Geometry3DParallelCameraMatrix::toString(std::string name =
                                                         "P") const;

} // namespace geometry
} // namespace KCT
