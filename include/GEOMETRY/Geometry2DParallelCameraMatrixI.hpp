#pragma once
// Logging
#include <plog/Log.h>

// Standard libraries
#include <array>
#include <iomanip>
#include <iostream>

// Internal libraries

namespace KCT {
namespace geometry {
    /**
     *Class to represent projection matrices of parallel rays with one detector coordinate.
     *
     * The projection matrix is a 1x4 row:
     *   PX = a0*x0 + a1*x1 + a2*x2 + a3
     */
    class Geometry2DParallelCameraMatrixI
    {
    public:
        virtual ~Geometry2DParallelCameraMatrixI() = default;

        /**
         * Vector VR is a unit vector in the direction of incoming rays.
         * Vector VX is a vector orthogonal to VR in the direction of positive increments of PX.
         * VX has a size corresponding to the real detector pixel size.
         */
        virtual void directionVectorVR(double* vector3) const = 0;
        virtual void directionVectorVX(double* vector3) const = 0;
        virtual std::array<double, 3> directionVectorVR() const = 0;
        virtual std::array<double, 3> directionVectorVX() const = 0;

        /**
         * Size of the detector pixel in world coordinates.
         *
         * @return |VX|
         */
        virtual double pixelSpacing() const = 0;

        /**
         * Get minimal vector (x,y,z) that gets projected to PX.
         * It corresponds to projector crossing origin.
         *
         * @param PX Detector coordinate
         *
         * @return Vector projecting to PX.
         */
        virtual void backprojectToPosition(const double PX, double* vector3, double z = 0.0) const
            = 0;
        virtual std::array<double, 3> backprojectToPosition(const double PX, double z = 0.0) const
            = 0;

        /**
         * Projection to the detector from the world coordinates.
         *
         * @param x0
         * @param x1
         * @param x2
         * @param PX
         */
        virtual void project(const double x0, const double x1, const double x2, double* PX) const
            = 0;
        virtual void project(const float x0, const float x1, const float x2, float* PX) const = 0;

        virtual void projectionMatrixAsVector4(double* vector4) const = 0;
        virtual std::array<double, 4> projectionMatrixAsVector4() const = 0;

        /**
         * Returns projection vector [PX, PY, offset] for given z coordinate.
         */
        virtual void projectionMatrixPXAsVector3(double* vector3, double z = 0.0) const = 0;
        virtual std::array<double, 3> projectionMatrixPXAsVector3(double z = 0.0) const = 0;
    };
} // namespace geometry
} // namespace KCT
