#include "GEOMETRY/Geometry3DParallel.hpp"

namespace KCT {
namespace geometry {
    /**
     *Class to represent parallel ray geometry by means of 2x4 homogeneous matrix
     *and detector tilt.
     */
    Geometry3DParallel::Geometry3DParallel(const Geometry3DParallelCameraMatrix& pcm,
                                           const double cosDetectorTilt)
        : parallelCameraMatrix(pcm)
        , cosDetectorTilt(cosDetectorTilt)
    {
    }
    // Analogous to ASTRA parallel3d_vec
    Geometry3DParallel::Geometry3DParallel(const std::array<double, 3> rayDirection,
                                           const std::array<double, 3> detectorOrigin,
                                           const std::array<double, 3> VX,
                                           const std::array<double, 3> VY)
        : parallelCameraMatrix(rayDirection, detectorOrigin, VX, VY)
    {
        // Compute cosDetectorTilt
        std::array<double, 3> normalToDetector = vectorProduct(VX, VY);
        normalToDetector = normalizeVector<3>(normalToDetector);
        std::array<double, 3> unitRayDirection = parallelCameraMatrix.directionVectorVR();
        cosDetectorTilt = std::abs(vectorDotProduct<3>(normalToDetector, unitRayDirection));
    }

    bool Geometry3DParallel::operator==(const Geometry3DParallel& rhs)
    {
        if(cosDetectorTilt == rhs.cosDetectorTilt
           && parallelCameraMatrix == rhs.parallelCameraMatrix)
        {
            return true;
        } else
        {
            return false;
        }
    }

    // Implements Geometry3DParallelI
    double Geometry3DParallel::pixelSkew() const { return parallelCameraMatrix.pixelSkew(); }
    double Geometry3DParallel::pixelArea() const { return parallelCameraMatrix.pixelArea(); }
    double Geometry3DParallel::detectorTilt() const { return cosDetectorTilt; }
    void Geometry3DParallel::detectorTilt(double* scalar) const { *scalar = cosDetectorTilt; }
    void Geometry3DParallel::directionVectorVR(double* vector3) const
    {
        parallelCameraMatrix.directionVectorVR(vector3);
    }
    void Geometry3DParallel::directionVectorVX(double* vector3) const
    {
        parallelCameraMatrix.directionVectorVX(vector3);
    }
    void Geometry3DParallel::directionVectorVY(double* vector3) const
    {
        parallelCameraMatrix.directionVectorVY(vector3);
    }
    /**
     *
     * @return Unit vector in the direction of incomming rays or its oposite.
     */
    std::array<double, 3> Geometry3DParallel::directionVectorVR() const
    {
        return parallelCameraMatrix.directionVectorVR();
    }
    std::array<double, 3> Geometry3DParallel::directionVectorVX() const
    {
        return parallelCameraMatrix.directionVectorVX();
    }
    std::array<double, 3> Geometry3DParallel::directionVectorVY() const
    {
        return parallelCameraMatrix.directionVectorVY();
    }
    void Geometry3DParallel::backprojectToPosition(const double PX,
                                                   const double PY,
                                                   double* vector3) const
    {
        parallelCameraMatrix.backprojectToPosition(PX, PY, vector3);
    }

    std::array<double, 3> Geometry3DParallel::backprojectToPosition(const double PX,
                                                                    const double PY) const
    {
        return parallelCameraMatrix.backprojectToPosition(PX, PY);
    }

    void Geometry3DParallel::project(
        const double x0, const double x1, const double x2, double* PX, double* PY) const
    {
        parallelCameraMatrix.project(x0, x1, x2, PX, PY);
    }
    void Geometry3DParallel::project(
        const float x0, const float x1, const float x2, float* PX, float* PY) const
    {
        parallelCameraMatrix.project(x0, x1, x2, PX, PY);
    }

    void Geometry3DParallel::projectionMatrixAsVector8(double* vector8) const
    {
        parallelCameraMatrix.projectionMatrixAsVector8(vector8);
    }

    std::array<double, 8> Geometry3DParallel::projectionMatrixAsVector8() const
    {
        return parallelCameraMatrix.projectionMatrixAsVector8();
    }

    void Geometry3DParallel::projectionMatrixPXAsVector4(double* vector4) const
    {
        parallelCameraMatrix.projectionMatrixPXAsVector4(vector4);
    }

    std::array<double, 4> Geometry3DParallel::projectionMatrixPXAsVector4() const
    {
        return parallelCameraMatrix.projectionMatrixPXAsVector4();
    }

    void Geometry3DParallel::projectionMatrixPXAsVector3(double* vector3, double z) const
    {
        parallelCameraMatrix.projectionMatrixPXAsVector3(vector3, z);
    }

    std::array<double, 3> Geometry3DParallel::projectionMatrixPXAsVector3(double z) const
    {
        return parallelCameraMatrix.projectionMatrixPXAsVector3(z);
    }

    // Static functions
    // Helper ASTRA parallel3d
    Geometry3DParallel
    Geometry3DParallel::initializeFromParameters(const double detector_spacing_x,
                                                 const double detector_spacing_y,
                                                 const uint32_t projection_size_x,
                                                 const uint32_t projection_size_y,
                                                 const double angle)
    {
        Geometry3DParallelCameraMatrix parallelCameraMatrix
            = Geometry3DParallelCameraMatrix::initializeFromParameters(
                detector_spacing_x, detector_spacing_y, projection_size_x, projection_size_y,
                angle);
        double cosDetectorTilt = 1.0;
        return Geometry3DParallel(parallelCameraMatrix, cosDetectorTilt);
    }

    // Helper ASTRA parallel3d_vec
    Geometry3DParallel
    Geometry3DParallel::initializeFromParameters(const uint32_t projection_size_x,
                                                 const uint32_t projection_size_y,
                                                 const std::array<double, 3> rayDirection,
                                                 const std::array<double, 3> detectorCenter,
                                                 const std::array<double, 3> VX,
                                                 const std::array<double, 3> VY)
    {
        double x_shift = -0.5 * double(projection_size_x);
        double y_shift = -0.5 * double(projection_size_y);
        std::array<double, 3> detector_x_shift = multiplyVectorByConstant<3>(VX, x_shift);
        std::array<double, 3> detector_y_shift = multiplyVectorByConstant<3>(VY, y_shift);
        std::array<double, 3> detectorOrigin = vectorSum<3>(detectorCenter, detector_x_shift);
        detectorOrigin = vectorSum<3>(detectorOrigin, detector_y_shift);
        return Geometry3DParallel(rayDirection, detectorOrigin, VX, VY);
    }

    // Helper equidistant angles
    Geometry3DParallel
    Geometry3DParallel::initializeFromParameters(const double detector_spacing_x,
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

} // namespace geometry
} // namespace KCT
