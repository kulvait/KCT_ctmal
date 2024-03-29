#pragma once
// Logging
#include <plog/Log.h>

// Internal libraries
#include "MATRIX/Matrix.hpp"

namespace KCT {
namespace matrix {
    class SquareMatrix : public Matrix
    {
    public:
        SquareMatrix(uint32_t n)
            : Matrix(n, n)
        {
        }

        SquareMatrix(const Matrix& squaremat)
            : Matrix(squaremat)
        {
            if(squaremat.dimm() != squaremat.dimn())
            {
                std::string msg = io::xprintf("Square matrix must have the same m and n dimensions "
                                              "but they differ m=%d, n=%d.",
                                              squaremat.dimm(), squaremat.dimn());
                LOGE << msg;
                throw new std::runtime_error(msg);
            }
        }

        SquareMatrix(uint32_t n, const double(&initArray))
            : Matrix(n, n, initArray)
        {
        }

        SquareMatrix(uint32_t n, std::initializer_list<double> A_copy)
            : Matrix(n, n, A_copy)
        {
        }

    private:
    };

} // namespace matrix
} // namespace KCT
