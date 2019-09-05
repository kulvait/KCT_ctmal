#pragma once
// Logging
#include <plog/Log.h>

// Standard libraries
#include <array>
#include <iomanip>
#include <iostream>

// Internal libraries
#include "MATRIX/LUDoolittleForm.hpp"
#include "MATRIX/Matrix.hpp"
#include "MATRIX/RQFactorization.hpp"
#include "MATRIX/SquareMatrix.hpp"

namespace CTL {
namespace matrix {
    /**
     *Class to represent projection matrices.
     */
    class ProjectionMatrix : public Matrix
    {
    private:
        void computeSourcePosition();

    public:
        /**Constructs new ProjectionMatrix that is inicialized by zeros.*/
        ProjectionMatrix()
            : Matrix(3, 4)
        {
            computeSourcePosition();
        }
        /**Constructor from the double array*/
        ProjectionMatrix(const double (&initArray)[3 * 4])
            : Matrix(3, 4, initArray)
        {
            computeSourcePosition();
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
            computeSourcePosition();
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
        /**Get normal to detector. The normal to the detector will points from the detector towards source.*/
        std::array<double, 3> normalToDetector() const;
        /**Get normalized vector from source to detector ending at (px,py)*/
        std::array<double, 3> projectedToPosition(double px, double py) const;
        /**
         *
         * @return Normalized vector that multiples of can be added to the source position to get
         * the normal to the detector in the x pixel direction.
         */
        std::array<double, 3> tangentToDetectorXDirection() const;
        /**
         *
         * @return Normalized vector that multiples of can be added to the source position to get
         * the normal to the detector in the y pixel direction.
         */
        std::array<double, 3> tangentToDetectorYDirection() const;
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

    private:
        /**
         * Compute normalized vector with respect to l2 norm of a 3D vector.
         *
         * @param v
         *
         * @return
         */
        inline std::array<double, 3> normalizeVector(std::array<double, 3> v) const;
        /**
         * Compute l2 norm of a 3D vector.
         *
         * @param v
         *
         * @return
         */
        std::array<double, 4> reorthogonalize(std::array<double, 4> v,
                                              std::array<double, 4> og) const;
        inline double vectorNorm(std::array<double, 3> v) const;
        bool sourceComputed = false;
        std::array<double, 3> source;
        const double zeroPrecisionTolerance = 1e-10;
    };
} // namespace matrix
} // namespace CTL
