#include "GEOMETRY/Geometry3DParallelCameraMatrix.hpp"

namespace KCT::geometry {

Geometry3DParallelCameraMatrix::Geometry3DParallelCameraMatrix(const double* vector8)
{
    std::copy(vector8, vector8 + 8, std::begin(this->projectionMatrixVector));
}

Geometry3DParallelCameraMatrix::Geometry3DParallelCameraMatrix(
    std::initializer_list<double> vector8_copy)
{
    if(vector8_copy.size() != 8)
    {
        KCTERR(io::xprintf(
            "The class must be initialized with proper sized list of 8 elements representing "
            "rows of 2x4 matrix but the size is %d.",
            vector8_copy.size()));
    }
    std::copy(vector8_copy.begin(), vector8_copy.end(), std::begin(this->projectionMatrixVector));
}
/**Constructor from the 2x4 Matrix.*/
Geometry3DParallelCameraMatrix::Geometry3DParallelCameraMatrix(const Matrix& pm) {}
// ASTRA like initialization
Geometry3DParallelCameraMatrix::Geometry3DParallelCameraMatrix(
    const std::array<double, 3> rayDirection,
    const std::array<double, 3> detectorOrigin,
    const std::array<double, 3> VX,
    const std::array<double, 3> VY)
{
    // See
    // https://kulvait.github.io/KCT_doc/posts/working-with-kct-cbct-5-parallel-beam-geometry.html
    // for explanation
    std::array<double, 3> VX0
        = orthogonalPartOfVectorWithRespectToSecondVector<3>(VX, rayDirection);
    std::array<double, 3> VY0
        = orthogonalPartOfVectorWithRespectToSecondVector<3>(VY, rayDirection);
    double VX0_lensquare = vectorDotProduct<3>(VX0, VX0);
    double VY0_lensquare = vectorDotProduct<3>(VY0, VY0);
    std::array<double, 3> a = multiplyVectorByConstant<3>(VX0, 1.0 / VX0_lensquare);
    std::array<double, 3> b = multiplyVectorByConstant<3>(VY0, 1.0 / VY0_lensquare);
    double px0 = -vectorDotProduct<3>(a, detectorOrigin);
    double py0 = -vectorDotProduct<3>(b, detectorOrigin);
    std::copy(std::begin(a), std::end(a), std::begin(projectionMatrixVector));
    projectionMatrixVector[3] = px0;
    std::copy(std::begin(b), std::end(b), &projectionMatrixVector[4]);
    projectionMatrixVector[7] = py0;
}

bool Geometry3DParallelCameraMatrix::operator==(const Geometry3DParallelCameraMatrix& rhs)
{
    for(int i = 0; i != projectionMatrixVector.size(); i++)
    {
        if(this->projectionMatrixVector[i] != rhs.projectionMatrixVector[i])
        {
            return false;
        }
    }
    return true;
}

// Implements Geometry3DParallelCameraMatrixI
double Geometry3DParallelCameraMatrix::pixelSkew() const
{
    // I have to normalize VX and VY and compute their scalar product
    std::array<double, 3> VXN = normalizeVector<3>(directionVectorVX());
    std::array<double, 3> VYN = normalizeVector<3>(directionVectorVY());
    // Here without absolute value as the vectors VX and VY direction matters
    return vectorDotProduct<3>(VXN, VYN);
}

double Geometry3DParallelCameraMatrix::pixelArea() const
{
    // I have to normalize VX and VY and compute their scalar product
    std::array<double, 3> VX = directionVectorVX();
    std::array<double, 3> VY = directionVectorVY();
    std::array<double, 3> AreaVector = vectorProduct(VX, VY);
    // Here without absolute value as the vectors VX and VY direction matters
    return vectorNorm<3>(AreaVector);
}

void Geometry3DParallelCameraMatrix::directionVectorVR(double* vector3) const
{
    std::array<double, 3> VXN = normalizeVector<3>(directionVectorVX());
    std::array<double, 3> VYN = normalizeVector<3>(directionVectorVY());
    std::array<double, 3> VR = vectorProduct(VXN, VYN);
    VR = normalizeVector<3>(VR);
    std::copy(std::begin(VR), std::begin(VR) + 3, vector3);
}

void Geometry3DParallelCameraMatrix::directionVectorVX(double* vector3) const
{
    std::array<double, 3> a;
    std::copy(std::begin(projectionMatrixVector), std::begin(projectionMatrixVector) + 3,
              std::begin(a));
    double normSquared = vectorDotProduct<3>(a, a);
    a = multiplyVectorByConstant(a, 1.0 / normSquared);
    std::copy(std::begin(a), std::begin(a) + 3, vector3);
}
void Geometry3DParallelCameraMatrix::directionVectorVY(double* vector3) const
{
    std::array<double, 3> b;
    std::copy(std::begin(projectionMatrixVector) + 4, std::begin(projectionMatrixVector) + 7,
              std::begin(b));
    double normSquared = vectorDotProduct<3>(b, b);
    b = multiplyVectorByConstant(b, 1.0 / normSquared);
    std::copy(std::begin(b), std::begin(b) + 3, vector3);
}
std::array<double, 3> Geometry3DParallelCameraMatrix::directionVectorVR() const
{
    std::array<double, 3> VR;
    directionVectorVR(std::begin(VR));
    return VR;
}
std::array<double, 3> Geometry3DParallelCameraMatrix::directionVectorVX() const
{
    std::array<double, 3> VX;
    directionVectorVX(std::begin(VX));
    return VX;
}
std::array<double, 3> Geometry3DParallelCameraMatrix::directionVectorVY() const
{
    std::array<double, 3> VY;
    directionVectorVY(std::begin(VY));
    return VY;
}
void Geometry3DParallelCameraMatrix::backprojectToPosition(const double PX,
                                                           const double PY,
                                                           double* vector3) const
{
    std::array<double, 3> x = backprojectToPosition(PX, PY);
    std::copy(std::begin(x), std::end(x), vector3);
}
std::array<double, 3> Geometry3DParallelCameraMatrix::backprojectToPosition(const double PX,
                                                                            const double PY) const
{
    // https://kulvait.github.io/KCT_doc/posts/working-with-kct-cbct-5-parallel-beam-geometry.html
    std::array<double, 3> a, b;
    std::copy(std::begin(projectionMatrixVector), std::begin(projectionMatrixVector) + 3,
              std::begin(a));
    std::copy(std::begin(projectionMatrixVector) + 4, std::begin(projectionMatrixVector) + 7,
              std::begin(b));
    double aa, ab, bb;
    aa = vectorDotProduct<3>(a, a);
    ab = vectorDotProduct<3>(a, b);
    bb = vectorDotProduct<3>(b, b);

    Matrix X = SquareMatrix(2, { aa, ab, ab, bb });
    LUDoolittleForm lu = LUDoolittleForm::LUDecomposeDoolittle(X, zeroPrecisionTolerance);
    SquareMatrix INV = lu.inverseMatrix();
    double alpha, beta;
    double PXO, PYO;
    PXO = PX - projectionMatrixVector[3];
    PYO = PY - projectionMatrixVector[7];
    alpha = INV(0, 0) * PXO + INV(0, 1) * PYO;
    beta = INV(1, 0) * PXO + INV(1, 1) * PYO;
    std::array<double, 3> x, x1, x2;
    x1 = multiplyVectorByConstant(a, alpha);
    x2 = multiplyVectorByConstant(b, beta);
    x = vectorSum(x1, x2);
    return x;
}
void Geometry3DParallelCameraMatrix::project(
    const double x0, const double x1, const double x2, double* PX, double* PY) const
{
    *PX = projectionMatrixVector[3] + projectionMatrixVector[0] * x0
        + projectionMatrixVector[1] * x1 + projectionMatrixVector[2] * x2;
    *PY = projectionMatrixVector[7] + projectionMatrixVector[4] * x0
        + projectionMatrixVector[5] * x1 + projectionMatrixVector[6] * x2;
}
void Geometry3DParallelCameraMatrix::project(
    const float x0, const float x1, const float x2, float* PX, float* PY) const
{
    *PX = projectionMatrixVector[3] + projectionMatrixVector[0] * x0
        + projectionMatrixVector[1] * x1 + projectionMatrixVector[2] * x2;
    *PY = projectionMatrixVector[7] + projectionMatrixVector[4] * x0
        + projectionMatrixVector[5] * x1 + projectionMatrixVector[6] * x2;
}

void Geometry3DParallelCameraMatrix::projectionMatrixAsVector8(double* vector8) const
{
    std::copy(std::begin(projectionMatrixVector), std::begin(projectionMatrixVector) + 8, vector8);
}
std::array<double, 8> Geometry3DParallelCameraMatrix::projectionMatrixAsVector8() const
{
    std::array<double, 8> pmv;
    projectionMatrixAsVector8(std::begin(pmv));
    return pmv;
}
std::string Geometry3DParallelCameraMatrix::toString(std::string name) const
{
    Matrix pmv(2, 4, std::begin(projectionMatrixVector));
    return pmv.toString(name);
}
// Static functions
// Helper ASTRA parallel3d
Geometry3DParallelCameraMatrix
Geometry3DParallelCameraMatrix::initializeFromParameters(const double detector_spacing_x,
                                                         const double detector_spacing_y,
                                                         const uint32_t projection_size_x,
                                                         const uint32_t projection_size_y,
                                                         const double omega)
{
    std::array<double, 8> pvm;
    std::array<double, 3> a, b, VR, VX, VY;
    // Astra toolbox consistent
    // See also
    // https://github.com/astra-toolbox/astra-toolbox/blob/fa2ec619edb2994e828897e80c06e7fb35c55c44/src/ParallelProjectionGeometry3D.cpp#L178
    VR[0] = std::cos(omega);
    VR[1] = std::sin(omega);
    VR[2] = 0.0;
    // Consistent with standard CBCT, see https://arxiv.org/pdf/2110.09841.pdf, Fig. 3
    VX[0] = std::sin(omega);
    VX[1] = -std::cos(omega);
    VX[2] = 0.0;
    a = multiplyVectorByConstant(VX, 1.0 / detector_spacing_x);
    VX = multiplyVectorByConstant(VX, detector_spacing_x);
    VY[0] = 0.0;
    VY[1] = 0.0;
    VY[2] = -detector_spacing_y;
    b[0] = 0.0;
    b[1] = 0.0;
    b[2] = -1.0 / detector_spacing_y;
    std::copy(std::begin(a), std::end(a), std::begin(pvm));
    pvm[3] = 0.5 * projection_size_x - 0.5;
    std::copy(std::begin(a), std::end(a), std::begin(pvm) + 4);
    pvm[7] = 0.5 * projection_size_y - 0.5;
    return Geometry3DParallelCameraMatrix(std::begin(pvm));
}

// Helper ASTRA parallel3d_vec
Geometry3DParallelCameraMatrix
Geometry3DParallelCameraMatrix::initializeFromParameters(const uint32_t projection_size_x,
                                                         const uint32_t projection_size_y,
                                                         const std::array<double, 3> rayDirection,
                                                         const std::array<double, 3> detectorCenter,
                                                         const std::array<double, 3> VX,
                                                         const std::array<double, 3> VY)
{
    double x_shift = -0.5 * (double(projection_size_x) - 1.0);
    double y_shift = -0.5 * (double(projection_size_y) - 1.0);
    std::array<double, 3> detector_x_shift = multiplyVectorByConstant<3>(VX, x_shift);
    std::array<double, 3> detector_y_shift = multiplyVectorByConstant<3>(VY, y_shift);
    std::array<double, 3> detectorOrigin = vectorSum<3>(detectorCenter, detector_x_shift);
    detectorOrigin = vectorSum<3>(detectorOrigin, detector_y_shift);
    return Geometry3DParallelCameraMatrix(rayDirection, detectorOrigin, VX, VY);
}

// Helper equidistant angles
Geometry3DParallelCameraMatrix
Geometry3DParallelCameraMatrix::initializeFromParameters(const double detector_spacing_x,
                                                         const double detector_spacing_y,
                                                         const uint32_t projection_size_x,
                                                         const uint32_t projection_size_y,
                                                         const uint32_t anglesCount,
                                                         const uint32_t angleNum)
{
    double angleIncrement = pi / anglesCount;
    double angle = angleNum * angleIncrement;
    return initializeFromParameters(detector_spacing_x, detector_spacing_y, projection_size_x,
                                    projection_size_y, angle);
}

} // namespace KCT::geometry
