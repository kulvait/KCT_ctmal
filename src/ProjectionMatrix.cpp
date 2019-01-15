// Logging on the top
#include "MATRIX/ProjectionMatrix.hpp"

namespace CTL {
namespace matrix {

    SquareMatrix<3> ProjectionMatrix::colSubMatrix(int j)
    {
        if(j < 0 || j > 3)
        {
            std::string errMsg = io::xprintf(
                "Row must be specified in the range [0,3], specified i=%d is out of range.", j);
            LOGE << errMsg;
            throw std::runtime_error(errMsg);
        }
        SquareMatrix<3> out;
        for(int i = 0; i != 3; i++)
        {
            for(int jtm = 0; jtm != 4; jtm++)
            {
                if(j > jtm)
                {
                    out(i, jtm) = (*this)(i, jtm);
                } else if(j == jtm)
                {
                    continue;
                } else
                {
                    out(i, jtm - 1) = (*this)(i, jtm);
                }
            }
        }
        return (out);
    }

    template <uint N>
    double determinant(const SquareMatrix<N>& M)
    {
        LUDoolittleForm<N> lu = LUDoolittleForm<N>::LUDecomposeDoolittle(M, 0.000001);
        return lu.getDeterminant();
    }

    std::array<double, 3> ProjectionMatrix::sourcePosition()
    {
        double divisor = -determinant<3>(this->colSubMatrix(3));
        std::array<double, 3> source;
        source[0] = determinant<3>(this->colSubMatrix(0)) / divisor;
        source[1] = -determinant<3>(this->colSubMatrix(1)) / divisor;
        source[2] = determinant<3>(this->colSubMatrix(2)) / divisor;
        return source;
    }

    /**
     *With a projection matrix provided by Siemens, I construct the vector that points to the 0 in
     *the world coordinates.
     */
    std::array<double, 3> ProjectionMatrix::normalToDetector()
    {
        std::array<double, 3> normal;
        normal[0] = (*this)(2, 0);
        normal[1] = (*this)(2, 1);
        normal[2] = (*this)(2, 2);
        double normvec
            = std::sqrt(normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]);
        normal[0] /= normvec;
        normal[1] /= normvec;
        normal[2] /= normvec;
        return normal;
    }

    /**
     *Actual projection
     */
    void ProjectionMatrix::project(double x, double y, double z, double* px, double* py)
    {
        double divideByMe
            = x * (*this)(2, 0) + y * (*this)(2, 1) + z * (*this)(2, 2) + (*this)(2, 3);
        *px = (x * (*this)(0, 0) + y * (*this)(0, 1) + z * (*this)(0, 2) + (*this)(0, 3))
            / divideByMe;
        *py = (x * (*this)(1, 0) + y * (*this)(1, 1) + z * (*this)(1, 2) + (*this)(1, 3))
            / divideByMe;
    }

    /**
     *Actual projection
     */
    void ProjectionMatrix::project(float x, float y, float z, float* px, float* py)
    {
        double divideByMe
            = x * (*this)(2, 0) + y * (*this)(2, 1) + z * (*this)(2, 2) + (*this)(2, 3);
        double pxd = (x * (*this)(0, 0) + y * (*this)(0, 1) + z * (*this)(0, 2) + (*this)(0, 3))
            / divideByMe;
        double pyd = (x * (*this)(1, 0) + y * (*this)(1, 1) + z * (*this)(1, 2) + (*this)(1, 3))
            / divideByMe;
        *px = float(pxd);
        *py = float(pyd);
    }

} // namespace util
} // namespace CTL
