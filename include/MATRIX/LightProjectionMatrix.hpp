#pragma once
// Logging
#include <plog/Log.h>

//External libs
#include <mkl.h>


//Internal libs
#include "MATRIX/utils.hpp"
#include "MATRIX/CameraI.hpp"
#include "MATRIX/ProjectionMatrix.hpp"

namespace KCT {
namespace matrix {
    /**
     *Class to represent projection matrices.
     */
    class LightProjectionMatrix : public CameraI
    {
    public:
        LightProjectionMatrix(const double* sourcePosition,const double* vector9);
        LightProjectionMatrix(const ProjectionMatrix& pm);
        // Implements CameraI
        std::array<double, 3> sourcePosition() const;
        double pixelSkew() const;
        std::array<double, 3> directionVectorVN() const;
        std::array<double, 3> directionVectorVX() const;
        std::array<double, 3> directionVectorVY() const;
        std::array<double, 2> principalRayProjection() const;
        std::array<double, 3> normalToDetector() const;
        std::array<double, 3> directionToPosition(const double pi, const double pj) const;
        void project(const double x0, const double x1, const double x2, double* pi, double* pj) const;
        void project(const float x0, const float x1, const float x2, float* pi, float* pj) const;
        void project0(const double X0, const double X1, const double X2, double* pi, double* pj) const;
        void project0( const float X0, const float X1, const float X2, float* pi, float* pj) const;
        std::array<double, 2> focalLength() const;
        std::array<double, 2> pixelSizes(const double sourceToDetector) const;
        double sourceToDetectorFromPX(const double PX) const;
        double sourceToDetectorFromPY(const double PY) const;
        void sourcePosition(double* vector3) const;
        void directionVectorVN(double* vector3) const;
        void directionVectorVX(double* vector3) const;
        void directionVectorVY(double* vector3) const;
        void normalToDetector(double* vector3) const;
        void principalRayProjection(double* vector2) const;
        void focalLength(double* vector2) const;
        void pixelSizes(const double sourceToDetector, double* vector2) const;
        void directionToPosition(const double pi, const double pj, double* vector3) const;
        void projectionMatrixAsVector12(double* vector12) const;
        void projectionMatrixAsVector9(double* vector9) const;
        void inverseProjectionMatrixAsVector16(double* vector16) const;
        void inverseProjectionMatrixAsVector9(double* vector9) const;
        std::array<double, 12> projectionMatrixAsVector12() const;
        std::array<double, 9> projectionMatrixAsVector9() const;
        std::array<double, 16> inverseProjectionMatrixAsVector16() const;
        std::array<double, 9> inverseProjectionMatrixAsVector9() const;

    private:
        std::array<double, 3> source;
        std::array<double, 9> matrix;//Projection matrix without fourth row
    };

} // namespace matrix
} // namespace KCT
