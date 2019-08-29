#pragma once
// Logging
#include <plog/Log.h>

// Standard libraries
#include <array>
#include <iomanip>
#include <iostream>

// Internal libraries
#include "MATRIX/Matrix.hpp"
#include "MATRIX/RQFactorization.hpp"
#include "MATRIX/SquareMatrix.hpp"
#include "MATRIX/LUDoolittleForm.hpp"

namespace CTL {
namespace matrix {
    /**
     *Class to represent projection matrices.
     */
    class ProjectionMatrix : public Matrix
    {
    public:
        /**Constructs new ProjectionMatrix that is inicialized by zeros.*/
        ProjectionMatrix()
            : Matrix(3, 4){};
        /**Constructor from the double array*/
        ProjectionMatrix(const double (&initArray)[3 * 4])
            : Matrix(3, 4, initArray)
        {
        }
        /**Constructor from the 3x4 Matrix.*/
        ProjectionMatrix(const Matrix& pm)
            : Matrix(pm)
        {
            if(pm.dimm() != 3 || pm.dimn() != 4)
            {
                std::string msg
                    = io::xprintf("Projection matrix must have dimension 3x4 and not %dx%d.",
                                  pm.dimm(), pm.dimn());
                LOGE << msg;
                throw new std::runtime_error(msg);
            }
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
        std::array<double, 3> sourcePosition() const;
        /**Get normal to detector ending at 0*/
        std::array<double, 3> normalToDetector();
        /*Get 3x3 submatrix of projection matrix, where i-th row is removed.*/
        SquareMatrix colSubMatrix(int i) const;
        /**Compute projection of volume point to the projector point*/
        void project(double x, double y, double z, double* px, double* py);
        void project(float x, float y, float z, float* px, float* py);

        /**Compute the projection matrix with origin shifted by (x,y) detector pixels
         *
         *The projection matrix is multiplied from the left by the matrix
         *1 0 -x
         *0 1 -y
         *0 0 1
         */
        ProjectionMatrix shiftDetectorOrigin(double x, double y) const;
        std::string toString(std::string name = "P") const;
    };
} // namespace matrix
} // namespace CTL
