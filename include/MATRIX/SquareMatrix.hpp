#pragma once
// Logging
#include <plog/Log.h>

// Internal libraries
#include "matrix.h"

namespace CTL {
namespace util {
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

} // namespace util
} // namespace CTL
