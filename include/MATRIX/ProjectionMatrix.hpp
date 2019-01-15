#pragma once
// Logging
#include <plog/Log.h>

// Standard libraries
#include <array>

// Internal libraries
#include "MATRIX/LUDoolittleForm.hpp"
#include "MATRIX/SquareMatrix.hpp"
#include "MATRIX/Matrix.hpp"

namespace CTL {
namespace matrix {
    /**
     *Class to represent projection matrices.
     */
    class ProjectionMatrix : public Matrix<3, 4>
    {
    public:
        /**Constructs new ProjectionMatrix that is inicialized by zeros.*/
        ProjectionMatrix()
            : Matrix<3, 4>(){};
        /**Constructor from the double array*/
        ProjectionMatrix(const double (&initArray)[3 * 4])
            : Matrix<3, 4>(initArray)
        {
        }
        /**Constructor from the 3x4 Matrix.*/
        ProjectionMatrix(const Matrix<3, 4>& pm)
            : Matrix<3, 4>(pm)
        {
        }
        /**Equality test*/

        bool operator==(const ProjectionMatrix& rhs)
        {
            for(int i = 0; i != 3; i++)
            {
                for(int j = 0; j != 4; j++)
                {
                    if((*this)(i, j) != (rhs)(i, j))
                        return false;
                }
            }
            return true;
        }
        /**Get source position*/
        std::array<double, 3> sourcePosition();
        /**Get normal to detector ending at 0*/
        std::array<double, 3> normalToDetector();
        /*Get 3x3 submatrix of projection matrix, where i-th row is removed.*/
        SquareMatrix<3> colSubMatrix(int i);
        /**Compute projection of volume point to the projector point*/
        void project(double x, double y, double z, double* px, double* py);
        void project(float x, float y, float z, float* px, float* py);
    };
} // namespace util
} // namespace CTL
