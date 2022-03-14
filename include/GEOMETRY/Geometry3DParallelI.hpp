#pragma once
// Logging
#include <plog/Log.h>

// Standard libraries
#include <array>
#include <iomanip>
#include <iostream>

// Internal libraries

namespace KCT {
namespace matrix {
    /**
     *Class to represent parallel ray geometry in 3D.
     */
    class Geometry3DParallelI
    {
    public:
        /**
         * @return Cos of the angle between VX and VY, dot product of unit vectors.
         */
        virtual double pixelSkew() const = 0;

        /**
         * @return Absolute value of cos of the angle between VR and normalToDetector, dot product of unit vectors.
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
         * Normal to the detector.
         *
         * @return Normal to detector pointing towards source.
         */
        virtual void normalToDetector(double* vector3) const = 0;
        virtual std::array<double, 3> normalToDetector() const = 0;

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
        backprojectToPosition(const double pi, const double pj, double* vector3) const = 0;
        virtual std::array<double, 3> backprojectToPosition(const double pi,
                                                          const double pj) const = 0;

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
            const double x0, const double x1, const double x2, double* pi, double* pj) const = 0;
        virtual void
        project(const float x0, const float x1, const float x2, float* pi, float* pj) const = 0;

        virtual void projectionMatrixAsVector8(double* vector8) const = 0;
        virtual std::array<double, 8> projectionMatrixAsVector8() const = 0;
    };
} // namespace matrix
} // namespace KCT
