#include "GEOMETRY/Geometry2DParallel.hpp"

namespace KCT {
namespace geometry {

    /**
     *Class to represent parallel ray geometry by means of 1x4 homogeneous matrix.
     */
    Geometry2DParallel::Geometry2DParallel(const Geometry2DParallelCameraMatrix& pcm)
        : parallelCameraMatrix(pcm)
    {
    }

    Geometry2DParallel::Geometry2DParallel(const std::array<double, 3>& rayDirection,
                                           const std::array<double, 3>& detectorOrigin,
                                           const std::array<double, 3>& VX)
        : parallelCameraMatrix(rayDirection, detectorOrigin, VX)
    {
    }

    bool Geometry2DParallel::operator==(const Geometry2DParallel& rhs) const
    {
        if(parallelCameraMatrix == rhs.parallelCameraMatrix)
        {
            return true;
        } else
        {
            return false;
        }
    }

    double Geometry2DParallel::pixelSpacing() const { return parallelCameraMatrix.pixelSpacing(); }

    void Geometry2DParallel::directionVectorVR(double* vector3) const
    {
        parallelCameraMatrix.directionVectorVR(vector3);
    }

    void Geometry2DParallel::directionVectorVX(double* vector3) const
    {
        parallelCameraMatrix.directionVectorVX(vector3);
    }

    std::array<double, 3> Geometry2DParallel::directionVectorVR() const
    {
        return parallelCameraMatrix.directionVectorVR();
    }

    std::array<double, 3> Geometry2DParallel::directionVectorVX() const
    {
        return parallelCameraMatrix.directionVectorVX();
    }

    void Geometry2DParallel::backprojectToPosition(const double PX, double* vector3) const
    {
        parallelCameraMatrix.backprojectToPosition(PX, vector3);
    }

    std::array<double, 3> Geometry2DParallel::backprojectToPosition(const double PX) const
    {
        return parallelCameraMatrix.backprojectToPosition(PX);
    }

    void Geometry2DParallel::project(const double x0, const double x1, const double x2,
                                     double* PX) const
    {
        parallelCameraMatrix.project(x0, x1, x2, PX);
    }

    void Geometry2DParallel::projectionMatrixAsVector4(double* vector4) const
    {
        parallelCameraMatrix.projectionMatrixAsVector4(vector4);
    }

    std::array<double, 4> Geometry2DParallel::projectionMatrixAsVector4() const
    {
        return parallelCameraMatrix.projectionMatrixAsVector4();
    }

    void Geometry2DParallel::projectionMatrixPXAsVector3(double* vector3, double z) const
    {
        parallelCameraMatrix.projectionMatrixPXAsVector3(vector3, z);
    }

    std::array<double, 3> Geometry2DParallel::projectionMatrixPXAsVector3(double z) const
    {
        return parallelCameraMatrix.projectionMatrixPXAsVector3(z);
    }

    // Static functions

    Geometry2DParallel
    Geometry2DParallel::initializeFromParameters(const double detector_spacing_x,
                                                 const uint32_t projection_size_x,
                                                 const double angle)
    {
        Geometry2DParallelCameraMatrix parallelCameraMatrix
            = Geometry2DParallelCameraMatrix::initializeFromParameters(detector_spacing_x,
                                                                       projection_size_x, angle);
        return Geometry2DParallel(parallelCameraMatrix);
    }

    Geometry2DParallel
    Geometry2DParallel::initializeFromParameters(const uint32_t projection_size_x,
                                                 const std::array<double, 3> rayDirection,
                                                 const std::array<double, 3> detectorCenter,
                                                 const std::array<double, 3> VX)
    {
        double x_shift = -0.5 * (double(projection_size_x) - 1.0);
        std::array<double, 3> detector_x_shift = multiplyVectorByConstant<3>(VX, x_shift);
        std::array<double, 3> detectorOrigin = vectorSum<3>(detectorCenter, detector_x_shift);
        return Geometry2DParallel(rayDirection, detectorOrigin, VX);
    }

    Geometry2DParallel
    Geometry2DParallel::initializeFromParameters(const double detector_spacing_x,
                                                 const uint32_t projection_size_x,
                                                 const uint32_t anglesCount,
                                                 const uint32_t angleNum)
    {
        double angleIncrement = pi / anglesCount;
        double angle = angleNum * angleIncrement;
        return initializeFromParameters(detector_spacing_x, projection_size_x, angle);
    }

} // namespace geometry
} // namespace KCT
