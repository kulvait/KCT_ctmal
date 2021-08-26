#pragma once
// Logging
#include <plog/Log.h>

// Standard libraries
#include <array>
#include <iomanip>
#include <iostream>

// Internal libraries
#include "MATRIX/CameraI.hpp"
#include "MATRIX/LUDoolittleForm.hpp"
#include "MATRIX/Matrix.hpp"
#include "MATRIX/RQFactorization.hpp"
#include "MATRIX/SquareMatrix.hpp"
#include "MATRIX/utils.hpp"

namespace KCT {
namespace matrix {
    /**
     *Class to represent projection matrices.
     */
    class ProjectionMatrix : public Matrix, CameraI
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
        // Implements CameraI
        std::array<double, 3> sourcePosition() const;
        double pixelSkew() const;
        std::array<double, 3> directionVectorVN() const;
        std::array<double, 3> directionVectorVX() const;
        std::array<double, 3> directionVectorVY() const;
        std::array<double, 2> principalRayProjection() const;
        std::array<double, 3> normalToDetector() const;
        std::array<double, 3> directionToPosition(const double pi, const double pj) const;
        void
        project(const double x0, const double x1, const double x2, double* pi, double* pj) const;
        void project(const float x0, const float x1, const float x2, float* pi, float* pj) const;
        void
        project0(const double X0, const double X1, const double X2, double* pi, double* pj) const;
        void project0(const float X0, const float X1, const float X2, float* pi, float* pj) const;
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
        bool sourceComputed = false;
        std::array<double, 3> source;
        const double zeroPrecisionTolerance = 1e-10;
    };
} // namespace matrix
} // namespace KCT
