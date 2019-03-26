#pragma once
// Logging, testing
#include <plog/Log.h>

// External libraries
#include <cmath>
#include <memory>

// Intel MKL
#include "mkl.h"

// Internal libraries
#include "MATRIX/Matrix.hpp"
#include "MATRIX/SquareMatrix.hpp"
#include "stringFormatter.h"

namespace CTL {
namespace matrix {
    /**RQ factorization is a Gramm Shmidt orthogonalization of rows of A from the last row.
     *
     *We require that m <= n since n is maximal dimension of m.
     *A = [R 0] [Q_1] = R Q_1, where R is mxm and Q_1 is mxn
     *          [Q_2]
     *We perform slim factorization to produce just R and Q_1.
     *See https://en.wikipedia.org/wiki/QR_decomposition
     */
    class RQFactorization
    {
    public:
        RQFactorization()
        {
            A = nullptr;
            R_factor = nullptr;
            tau = nullptr;
        }

        void factorize(std::shared_ptr<Matrix> A);
        std::shared_ptr<Matrix> getRMatrix(); // m*m
        std::shared_ptr<Matrix> getQMatrix(); // m*n

    private:
        std::shared_ptr<Matrix> A = nullptr; // m*n
        std::shared_ptr<Matrix> R_factor = nullptr; // Upper diagonal result matrix m*n
        std::shared_ptr<Matrix> tau
            = nullptr; // Elementary reflections, see
                       // https://software.intel.com/en-us/mkl-developer-reference-c-orgrq m*1
    };
} // namespace matrix
} // namespace CTL
