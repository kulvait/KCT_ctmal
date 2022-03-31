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
     *Class to represent parallel ray geometry in 3D.
     */
    class Geometry3DParallelI
    {
    public:
        /**
         * Cosine of the angle between PX and PY direction on the detector in word coordinates, it
         * is 0.0 for rectangular pixels.
         *
         * @return Cos of the angle between VX and VY, dot product of unit vectors, 0.0 for
         * rectangular pixels.
         */
        virtual double pixelSkew() const = 0;

        /**
         * Area of the pixel on the detector.
         *
         * @return |PX x PY| or |PX| |PY| sin (skew)
         */
        virtual double pixelArea() const = 0;

        /**
         * Cosine of the angle between the detector and surface orthogonal to incoming rays, for
         * detector orthogonal to the rays it is 1.0.
         *
         * @return Absolute value of cosine of the normal vectors, often 1.0
         */
        virtual double detectorTilt() const = 0;
        virtual void detectorTilt(double* scalar) const = 0;

        /**
         * Vector VR is a unit vector in the direction of incomming rays.
         * Vector VX is a vector orthogonal to VN in the direction of positive increments of PX.
         * Vector VY is a vector orthogonal to VN in the direction of positive increments of PY.
         * For rectangular voxels VX and VY are orthogonal.
         * VX and VY have a size that is corresponding to the real pixel size.
         */
        virtual void directionVectorVR(double* vector3) const = 0;
        virtual void directionVectorVX(double* vector3) const = 0;
        virtual void directionVectorVY(double* vector3) const = 0;
        virtual std::array<double, 3> directionVectorVR() const = 0;
        virtual std::array<double, 3> directionVectorVX() const = 0;
        virtual std::array<double, 3> directionVectorVY() const = 0;

        /**
         * Get minimal vector (x,y,z) that gets projected to (PX, PY).
         * It corresponds to projector crossing origin.
         *
         * @param pi
         * @param pj
         *
         * @return Vector projecting to (pi,pj).
         */
        virtual void
        backprojectToPosition(const double PX, const double PY, double* vector3) const = 0;
        virtual std::array<double, 3> backprojectToPosition(const double PX,
                                                            const double PY) const = 0;

        /**
         * Projection to the detector from the word coordinates.
         *
         * @param x
         * @param y
         * @param z
         * @param pi
         * @param pj
         */
        virtual void project(
            const double x0, const double x1, const double x2, double* PX, double* PY) const = 0;
        virtual void
        project(const float x0, const float x1, const float x2, float* PX, float* PY) const = 0;

        virtual void projectionMatrixAsVector8(double* vector8) const = 0;
        virtual std::array<double, 8> projectionMatrixAsVector8() const = 0;
    };
} // namespace geometry
} // namespace KCT
