#pragma once
// Logging
#include <plog/Log.h>

// Standard libraries
#include <array>
#include <iomanip>
#include <iostream>

// Internal libraries

namespace CTL {
namespace matrix {
    /**
     *Class to represent camera and projective system.
     */
    class CameraI
    {
    public:
        virtual void sourcePosition(double* vector3) const = 0;
        /**
         * @return Source position in a word coordinates.
         */
        virtual std::array<double, 3> sourcePosition() const = 0;

        /**
         * Let's move the word coordinates (x,y,z) to the source related coordinates (X, Y, Z) =
         * (x-sx, y-sy, z-sz) Vector VN is a unit vector orthogonal to the detector such that (VNX,
         * VNY, VNZ).(-sx, -sy, -sz) >= 0. So when we assume that the volume is oriented around
         * (x,y,z)=(0,0,0) that is the directional vector of the principal ray towards the volume.
         * Vector VX is a vector orthogonal to VN in the direction of positive increments of PX.
         * Vector VY is a vector orthogonal to VN in the direction of positive increments of PY.
         * For rectangular voxels VX and VY are orthogonal.
         *
         * @return Cos of the angle between VX and VY, dot product of unit vectors.
         */
        virtual double pixelSkew() const = 0;

        /**
         * Let's move the word coordinates (x,y,z) to the source related coordinates (X, Y, Z) =
         * (x-sx, y-sy, z-sz) Vector VN is a unit vector orthogonal to the detector such that (VNX,
         * VNY, VNZ).(-sx, -sy, -sz) >= 0. So when we assume that the volume is oriented around
         * (x,y,z)=(0,0,0) that is the directional vector of the principal ray towards the volume.
         * Vector VX is a vector orthogonal to VN in the direction of positive increments of PX.
         * Vector VY is a vector orthogonal to VN in the direction of positive increments of PY.
         * For rectangular voxels VX and VY are orthogonal.
         * VX and VY have a size that is focal length.
         */
        virtual void directionVectorVN(double* vector3) const = 0;
        virtual void directionVectorVX(double* vector3) const = 0;
        virtual void directionVectorVY(double* vector3) const = 0;
        virtual std::array<double, 3> directionVectorVN() const = 0;
        virtual std::array<double, 3> directionVectorVX() const = 0;
        virtual std::array<double, 3> directionVectorVY() const = 0;

        /**
         * Normal to the detector that points from detector surface towards source
         *
         * Let's move the word coordinates (x,y,z) to the source related coordinates (X, Y, Z) =
         * (x-sx, y-sy, z-sz) Vector VN is a unit vector orthogonal to the detector such that (VNX,
         * VNY, VNZ)*(-sx, -sy, -sz) >= 0. So when we assume that the volume is oriented around
         * (x,y,z)=(0,0,0) that is the directional vector of the principal ray towards the volume.
         * Then normal to the detector will be -VN.
         *
         * @return Normal to detector pointing towards source.
         */
        virtual void normalToDetector(double* vector3) const = 0;
        virtual std::array<double, 3> normalToDetector() const = 0;

        /**
         * @return (pi, pj) is the position of the tangent line from source to the detector on the
         * detector.
         */
        virtual void principalRayProjection(double* vector2) const = 0;
        virtual std::array<double, 2> principalRayProjection() const = 0;

        /**
         * @return (fi, fj) that is focal length .
         */
        virtual void focalLength(double* vector2) const = 0;
        virtual std::array<double, 2> focalLength() const = 0;

        /**
         * @return Pixel spacings when sourceToDetector norm is given.
         */
        virtual void pixelSizes(const double sourceToDetector, double* vector2) const = 0;
        virtual std::array<double, 2> pixelSizes(const double sourceToDetector) const = 0;

        /**
         * @return Source to detector distance when px is given pixel spacing in x direction.
         */
        virtual double sourceToDetectorFromPX(const double PX) const = 0;

        /**
         * @return Source to detector distance when py is given pixel spacing in y direction.
         */
        virtual double sourceToDetectorFromPY(const double PY) const = 0;

        /**
         * Get normalized directional vector (X,Y,Z), such that (VNX, VNY, VNZ).(X,Y,Z) > 0 and that
         * s+X gets projected to (PX, PY)
         *
         * @param pi
         * @param pj
         *
         * @return Normalized vector from source to detector ending at (pi,pj).
         */
        virtual void
        directionToPosition(const double pi, const double pj, double* vector3) const = 0;
        virtual std::array<double, 3> directionToPosition(const double pi,
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
        /**
         * Projection to the detector from the source coordinates.
         *
         * @param X = x-sx
         * @param Y = y-sy
         * @param Z = z -sz
         * @param pi
         * @param pj
         */
        virtual void project0(
            const double X0, const double X1, const double X2, double* pi, double* pj) const = 0;
        virtual void
        project0(const float X0, const float X1, const float X2, float* pi, float* pj) const = 0;

        virtual void projectionMatrixAsVector12(double* vector12) const = 0;
        virtual std::array<double, 12> projectionMatrixAsVector12() const = 0;

        virtual void projectionMatrixAsVector9(double* vector9) const = 0;
        virtual std::array<double, 9> projectionMatrixAsVector9() const = 0;

        virtual void inverseProjectionMatrixAsVector16(double* vector16) const = 0;
        virtual std::array<double, 16>
        inverseProjectionMatrixAsVector16() const = 0;

        virtual void inverseProjectionMatrixAsVector9(double* vector9) const = 0;
        virtual std::array<double, 9> inverseProjectionMatrixAsVector9() const = 0;
    };
} // namespace matrix
} // namespace CTL
