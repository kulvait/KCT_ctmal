#pragma once
// Logging
#include <plog/Log.h>

// Standard libraries
#include <array>
#include <iomanip>
#include <iostream>

// Internal libraries
#include "GEOMETRY/Geometry3DParallelCameraMatrixI.hpp"

namespace KCT {
namespace geometry {
    /**
     *Class to represent projection matrices of parallel rays that are 2x4 matrices.
     */
    class Geometry3DParallelCameraMatrix : public Geometry3DParallelCameraMatrixI
    {
    private:
        void computeRayDirection();

    public:
        /**Constructs new ProjectionMatrix that is inicialized by zeros.*/
        ProjectionMatrix()
            : Matrix(2, 4)
        {
            computeRayDirection();
        }
        /**Constructor from the double array*/
        ProjectionMatrix(const double (&initArray)[2 * 4])
            : Matrix(2, 4, initArray)
        {
            computeRayDirection();
        }
        /**Constructor from the 2x4 Matrix.*/
        ProjectionMatrix(const Matrix& pm)
            : Matrix(pm)
        {
            if(pm.dimm() != 2 || pm.dimn() != 4)
            {
                std::string msg = io::xprintf(
                    "Parallel projection matrix must have dimension 2x4 and not %dx%d.", pm.dimm(),
                    pm.dimn());
                LOGE << msg;
                throw new std::runtime_error(msg);
            }
            computeRayDirection();
        }
        /**Equality test*/

        bool operator==(const Geometry3DParallelCameraMatrix& rhs)
        {
            for(int i = 0; i != 2; i++)
            {
                for(int j = 0; j != 4; j++)
                {
                    if((*this)(i, j) != (rhs)(i, j))
                        return false;
                }
            }
            return true;
        }
        // Implements Geometry3DParallelCameraMatrixI
        double pixelSkew() const;
        void directionVectorVR(double* vector3) const;
        void directionVectorVX(double* vector3) const;
        void directionVectorVY(double* vector3) const;
        std::array<double, 3> directionVectorVR() const;
        std::array<double, 3> directionVectorVX() const;
        std::array<double, 3> directionVectorVY() const;
        void
        backprojectToPosition(const double pi, const double pj, double* vector3) const;
        std::array<double, 3> backprojectToPosition(const double pi,
                                                            const double pj) const;
        void project(
            const double x0, const double x1, const double x2, double* pi, double* pj) const;
        void
        project(const float x0, const float x1, const float x2, float* pi, float* pj) const;

        void projectionMatrixAsVector8(double* vector8) const;
        std::array<double, 8> projectionMatrixAsVector8() const;
        std::string toString(std::string name = "P") const;

    private:
        std::array<double, 8> matrix; // Projection matrix
        const double zeroPrecisionTolerance = 1e-10;
    };
} // namespace geometry
} // namespace KCT
