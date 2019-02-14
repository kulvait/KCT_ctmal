#pragma once
// Logging
#include <plog/Log.h>

// Internal libraries
#include "MATRIX/Matrix.hpp"

namespace CTL {
namespace matrix {
    template <uint N>
    class SquareMatrix : public Matrix<N, N>
    {
    public:
        SquareMatrix()
            : Matrix<N, N>()
        {
        }

        SquareMatrix(const Matrix<N, N>& squaremat)
            : Matrix<N, N>(squaremat)
        {
        }

        SquareMatrix(const double (&initArray)[N * N])
            : Matrix<N, N>(initArray)
        {
        }

    private:
    };

} // namespace matrix
} // namespace CTL
