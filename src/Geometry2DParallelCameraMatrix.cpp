#include "GEOMETRY/Geometry2DParallelCameraMatrix.hpp"

namespace KCT::geometry {

Geometry2DParallelCameraMatrix::Geometry2DParallelCameraMatrix(const double* vector4)
{
    std::copy(vector4, vector4 + 4, std::begin(this->projectionMatrixVector));
}

Geometry2DParallelCameraMatrix::Geometry2DParallelCameraMatrix(
    std::initializer_list<double> vector4_copy)
{
    if(vector4_copy.size() != 4)
    {
        KCTERR(io::xprintf(
            "The class must be initialized with proper sized list of 4 elements representing "
            "rows of 1x4 matrix but the size is %d.",
            vector4_copy.size()));
    }
    std::copy(vector4_copy.begin(), vector4_copy.end(), std::begin(this->projectionMatrixVector));
}

// ASTRA-like initialization from geometry description
Geometry2DParallelCameraMatrix::Geometry2DParallelCameraMatrix(
    const std::array<double, 3> rayDirection,
    const std::array<double, 3> detectorOrigin,
    const std::array<double, 3> VX)
{
    std::array<double, 3> VX0
        = orthogonalPartOfVectorWithRespectToSecondVector<3>(VX, rayDirection);
    double VX0_lensquare = vectorDotProduct<3>(VX0, VX0);
    std::array<double, 3> a = multiplyVectorByConstant<3>(VX0, 1.0 / VX0_lensquare);
    double px0 = -vectorDotProduct<3>(a, detectorOrigin);

    std::copy(std::begin(a), std::end(a), std::begin(projectionMatrixVector));
    projectionMatrixVector[3] = px0;
}

bool Geometry2DParallelCameraMatrix::operator==(const Geometry2DParallelCameraMatrix& rhs) const
{
    for(uint32_t i = 0; i != projectionMatrixVector.size(); i++)
    {
        if(this->projectionMatrixVector[i] != rhs.projectionMatrixVector[i])
        {
            return false;
        }
    }
    return true;
}

// We use first two elements of projection matrix, third entry is in this geontery only to adjust
// offset of projection and is not used to determine direction of detector or rays.
void Geometry2DParallelCameraMatrix::directionVectorVX(double* vector3) const
{
    std::array<double, 2> a;
    std::copy(std::begin(projectionMatrixVector), std::begin(projectionMatrixVector) + 2,
              a.begin());
    std::array<double, 3> VX;
    double normSquared = vectorDotProduct<2>(a, a);
    if(normSquared < zeroPrecisionTolerance)
    {
        KCTERR("Can not determine detector direction VX from projection matrix.");
    }

    // Get vector in the same direction as a but with length equal to pixel spacing (i.e. norm of
    // VX)
    VX[0] = a[0] / normSquared;
    VX[1] = a[1] / normSquared;
    VX[2] = 0.0;

    std::copy(std::begin(VX), std::end(VX), vector3);
}

std::array<double, 3> Geometry2DParallelCameraMatrix::directionVectorVX() const
{
    std::array<double, 3> VX;
    directionVectorVX(std::begin(VX));
    return VX;
}

void Geometry2DParallelCameraMatrix::directionVectorVR(double* vector3) const
{
    std::array<double, 3> VXN = normalizeVector<3>(directionVectorVX());
    std::array<double, 3> VYN = { 0.0, 0.0, 1.0 }; // Normal to XY plane
    std::array<double, 3> VR = vectorProduct(VXN, VYN);
    VR = normalizeVector<3>(VR);
    std::copy(std::begin(VR), std::end(VR), vector3);
}

std::array<double, 3> Geometry2DParallelCameraMatrix::directionVectorVR() const
{
    std::array<double, 3> VR;
    directionVectorVR(std::begin(VR));
    return VR;
}

double Geometry2DParallelCameraMatrix::pixelSpacing() const
{
    std::array<double, 3> VX = directionVectorVX();
    return vectorNorm<3>(VX);
}

void Geometry2DParallelCameraMatrix::backprojectToPosition(const double PX,
                                                           double* vector3,
                                                           double z) const
{
    std::array<double, 3> x = backprojectToPosition(PX);
    std::copy(std::begin(x), std::end(x), vector3);
}

// We are looking to the position on the detector alpha a + PX0 = PX for given offset parametrized by z
std::array<double, 3> Geometry2DParallelCameraMatrix::backprojectToPosition(const double PX,
                                                                            double z) const
{
    std::array<double, 2> a;
    std::copy(std::begin(projectionMatrixVector), std::begin(projectionMatrixVector) + 2,
              a.begin());
    double aa = vectorDotProduct<2>(a, a);

    double OFFSETX = projectionMatrixVector[3] + projectionMatrixVector[2] * z;
    double PXO = PX - OFFSETX;
    if(aa < zeroPrecisionTolerance)
    {
        KCTERR("Can not backproject from singular projection matrix.");
    }
    double alpha = PXO / aa;
    std::array<double, 2> detectorPosition = multiplyVectorByConstant<2>(a, alpha);
    std::array<double, 3> position = { detectorPosition[0], detectorPosition[1], z };
    return position;
}

void Geometry2DParallelCameraMatrix::project(const double x0,
                                             const double x1,
                                             const double x2,
                                             double* PX) const
{
    *PX = projectionMatrixVector[3] + projectionMatrixVector[0] * x0
        + projectionMatrixVector[1] * x1 + projectionMatrixVector[2] * x2;
}

void Geometry2DParallelCameraMatrix::project(const float x0,
                                             const float x1,
                                             const float x2,
                                             float* PX) const
{
    *PX = projectionMatrixVector[3] + projectionMatrixVector[0] * x0
        + projectionMatrixVector[1] * x1 + projectionMatrixVector[2] * x2;
}

void Geometry2DParallelCameraMatrix::projectionMatrixAsVector4(double* vector4) const
{
    std::copy(std::begin(projectionMatrixVector), std::end(projectionMatrixVector), vector4);
}

std::array<double, 4> Geometry2DParallelCameraMatrix::projectionMatrixAsVector4() const
{
    std::array<double, 4> pmv;
    projectionMatrixAsVector4(std::begin(pmv));
    return pmv;
}

void Geometry2DParallelCameraMatrix::projectionMatrixPXAsVector3(double* vector3, double z) const
{
    std::copy(std::begin(projectionMatrixVector), std::begin(projectionMatrixVector) + 2, vector3);
    vector3[2] = projectionMatrixVector[3] + projectionMatrixVector[2] * z;
}

std::array<double, 3> Geometry2DParallelCameraMatrix::projectionMatrixPXAsVector3(double z) const
{
    std::array<double, 3> pmv;
    projectionMatrixPXAsVector3(std::begin(pmv), z);
    return pmv;
}

std::string Geometry2DParallelCameraMatrix::toString(std::string name) const
{
    Matrix pmv(1, 4, std::begin(projectionMatrixVector));
    return pmv.toString(name);
}

// Static functions

Geometry2DParallelCameraMatrix Geometry2DParallelCameraMatrix::initializeFromParameters(
    const double detector_spacing_x, const uint32_t projection_size_x, const double omega)
{
    std::array<double, 4> pvm;
    std::array<double, 3> a, VR, VX;

    // Analogous to 3D parallel in XY plane
    VR[0] = std::cos(omega);
    VR[1] = std::sin(omega);
    VR[2] = 0.0;

    VX[0] = std::sin(omega);
    VX[1] = -std::cos(omega);
    VX[2] = 0.0;

    a = multiplyVectorByConstant(VX, 1.0 / detector_spacing_x);

    std::copy(std::begin(a), std::end(a), std::begin(pvm));
    pvm[3] = 0.5 * projection_size_x - 0.5;

    return Geometry2DParallelCameraMatrix(std::begin(pvm));
}

Geometry2DParallelCameraMatrix
Geometry2DParallelCameraMatrix::initializeFromParameters(const uint32_t projection_size_x,
                                                         const std::array<double, 3> rayDirection,
                                                         const std::array<double, 3> detectorCenter,
                                                         const std::array<double, 3> VX)
{
    double x_shift = -0.5 * (double(projection_size_x) - 1.0);
    std::array<double, 3> detector_x_shift = multiplyVectorByConstant<3>(VX, x_shift);
    std::array<double, 3> detectorOrigin = vectorSum<3>(detectorCenter, detector_x_shift);
    return Geometry2DParallelCameraMatrix(rayDirection, detectorOrigin, VX);
}

Geometry2DParallelCameraMatrix
Geometry2DParallelCameraMatrix::initializeFromParameters(const double detector_spacing_x,
                                                         const uint32_t projection_size_x,
                                                         const uint32_t anglesCount,
                                                         const uint32_t angleNum)
{
    double angleIncrement = pi / anglesCount;
    double angle = angleNum * angleIncrement;
    return initializeFromParameters(detector_spacing_x, projection_size_x, angle);
}

} // namespace KCT::geometry
