#pragma once
// Logging, testing
#include <plog/Log.h>

// External libraries
#include <cmath>
#include <memory>

// Internal libraries
#include "MATRIX/Matrix.hpp"
#include "MATRIX/SquareMatrix.hpp"
#include "stringFormatter.h"

namespace CTL {
namespace matrix {
    class LUDoolittleForm
    {
    public:
        LUDoolittleForm(std::shared_ptr<SquareMatrix> LU,
                        std::shared_ptr<std::vector<uint32_t>> P,
                        bool oddSwapParity)
            : LU(LU)
            , P(P)
            , oddSwapParity(oddSwapParity)
        {
		n=LU->dimn();
        }

        double getDeterminant();
        double getLog10Determinant();
        SquareMatrix inverseMatrix();
        Matrix solve(Matrix vector); // Nx1 to Nx1
        SquareMatrix getLMatrix();
        SquareMatrix getUMatrix();
        SquareMatrix getPAMatrix();
        std::vector<uint32_t> getPermutation();
        bool getOddSwapParity();
        static LUDoolittleForm
        LUDecomposeDoolittle(const SquareMatrix& A, double minimalPivotSize);
        Matrix backwardSubstitute(Matrix b); // Nx1 to Nx1
    private:
        uint32_t n;
        std::shared_ptr<SquareMatrix> LU;
        std::shared_ptr<std::vector<uint32_t>> P;
        Matrix forwardSubstitute(Matrix b);
        static uint32_t findPivot(const SquareMatrix& A, const std::vector<uint32_t>& P, uint32_t k, double minimalPivotSize);
        bool oddSwapParity;
    };
} // namespace matrix
} // namespace CTL
