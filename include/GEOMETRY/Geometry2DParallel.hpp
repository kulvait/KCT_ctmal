#pragma once
// Logging
#include <plog/Log.h>

// Standard libraries
#include <array>
#include <cstdint>

// Internal libs
#include "GEOMETRY/Geometry2DParallelCameraMatrix.hpp"
#include "MATRIX/utils.hpp"

namespace KCT {
namespace geometry {

    class Geometry2DParallel
    {
    public:
        Geometry2DParallel(const Geometry2DParallelCameraMatrix& pcm);

        Geometry2DParallel(const std::array<double, 3>& rayDirection,
                           const std::array<double, 3>& detectorOrigin,
                           const std::array<double, 3>& VX);

        Geometry2DParallel(const std::array<double, 4>& projectionMatrix,
                           const std::array<double, 3>& rayDirection);

        bool operator==(const Geometry2DParallel& rhs) const;

        void directionVectorVR(double* vector3) const;
        void directionVectorVX(double* vector3) const;
        std::array<double, 3> directionVectorVR() const;
        std::array<double, 3> directionVectorVX() const;

        double pixelSpacing() const;

        void project(const double x0, const double x1, const double x2, double* PX) const;
        double project(const std::array<double, 3>& x) const;
        
        void projectionMatrixAsVector4(double* vector4) const;
        std::array<double, 4> projectionMatrixAsVector4() const;

        void projectionMatrixPXAsVector3(double* vector3, double z = 0.0) const;
        std::array<double, 3> projectionMatrixPXAsVector3(double z = 0.0) const;

        // Returns one point on the line corresponding to detector coordinate PX
        void backprojectToPosition(const double PX, double* vector3) const;
        std::array<double, 3> backprojectToPosition(const double PX) const;

        // Static helpers analogous to the 3D API
        static Geometry2DParallel initializeFromParameters(const double detector_spacing_x,
                                                           const uint32_t projection_size_x,
                                                           const double angle);

        static Geometry2DParallel
        initializeFromParameters(const uint32_t projection_size_x,
                                 const std::array<double, 3> rayDirection,
                                 const std::array<double, 3> detectorCenter,
                                 const std::array<double, 3> VX);

        static Geometry2DParallel initializeFromParameters(const double detector_spacing_x,
                                                           const uint32_t projection_size_x,
                                                           const uint32_t anglesCount,
                                                           const uint32_t angleNum);

    private:
        std::array<double, 4> projectionMatrixVector;
        std::array<double, 3> rayDirectionVector;
        static constexpr double zeroPrecisionTolerance = 1e-10;
        static constexpr double pi = 3.141592653589793238462643383279502884;
        Geometry2DParallelCameraMatrix parallelCameraMatrix;
    };

} // namespace geometry
} // namespace KCT
