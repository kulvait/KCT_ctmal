#pragma once
// Logging
#include <plog/Log.h>

// Standard libraries
#include <array>
#include <iomanip>
#include <iostream>

// Internal libraries
#include "MATRIX/Matrix.hpp"

// Interfaces
#include "GEOMETRY/Geometry3DParallelCameraMatrixI.hpp"

namespace KCT {
namespace geometry {

    using namespace KCT::matrix;

    /**
     *Class to represent projection matrices of parallel rays that are 2x4 matrices.
     */
    class Geometry3DParallelCameraMatrix : public Geometry3DParallelCameraMatrixI
    {
    private:
        void computeRayDirection();

    public:
        Geometry3DParallelCameraMatrix(const double* vector8);
        Geometry3DParallelCameraMatrix(const double (&initArray)[2 * 4]);
        /**Constructor from the 2x4 Matrix.*/
        Geometry3DParallelCameraMatrix(const Matrix& pm);
        // ASTRA like initialization
        Geometry3DParallelCameraMatrix(const std::array<double, 3> rayDirection,
                                       const std::array<double, 3> detectorOrigin,
                                       const std::array<double, 3> VX,
                                       const std::array<double, 3> VY);
        Geometry3DParallelCameraMatrix(const uint32_t anglesCount,
                                       const uint32_t angleNum,
                                       const double x_spacing,
                                       const double y_spacing);
        bool operator==(const Geometry3DParallelCameraMatrix& rhs)
        {
            for(int i = 0; i != 8; i++)
            {
                if(this->matrix[i] != rhs.matrix[i])
                    return false;
            }
            return true;
        }
        // Implements Geometry3DParallelCameraMatrixI
        double pixelSkew() const;

        void directionVectorVR(double* vector3) const;
        void directionVectorVX(double* vector3) const;
        void directionVectorVY(double* vector3) const;
        /**
         *
         * @return Unit vector in the direction of incomming rays or its oposite.
         */
        std::array<double, 3> directionVectorVR() const;
        std::array<double, 3> directionVectorVX() const;
        std::array<double, 3> directionVectorVY() const;
        void backprojectToPosition(const double pi, const double pj, double* vector3) const;
        std::array<double, 3> backprojectToPosition(const double pi, const double pj) const;
        void
        project(const double x0, const double x1, const double x2, double* pi, double* pj) const;
        void project(const float x0, const float x1, const float x2, float* pi, float* pj) const;

        void projectionMatrixAsVector8(double* vector8) const;
        std::array<double, 8> projectionMatrixAsVector8() const;
        std::string toString(std::string name = "P") const;

    private:
        std::array<double, 8> projectionMatrixVector; // Projection matrix
        const double zeroPrecisionTolerance = 1e-10;
    };
} // namespace geometry
} // namespace KCT
