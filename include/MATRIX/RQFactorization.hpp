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
    template <uint m, uint n>
    class RQFactorization
    {
    public:
        RQFactorization()
        {
            if(m > n)
            {
                io::throwerr("Initial matrix must be in the form m <= n but m=%d, n=%d.", m, n);
            }
            A = nullptr;
            R_factor = nullptr;
            tau = nullptr;
        }

        void factorize(std::shared_ptr<Matrix<m, n>> A);
        std::shared_ptr<Matrix<m, m>> getRMatrix();
        std::shared_ptr<Matrix<m, n>> getQMatrix();

    private:
        std::shared_ptr<Matrix<m, n>> A = nullptr;
        std::shared_ptr<Matrix<m, n>> R_factor = nullptr; // Upper diagonal result matrix
        std::shared_ptr<Matrix<m, 1>> tau
            = nullptr; // Elementary reflections, see
                       // https://software.intel.com/en-us/mkl-developer-reference-c-orgrq
    };

    template <uint m, uint n>
    void RQFactorization<m, n>::factorize(std::shared_ptr<Matrix<m, n>> A)
    {
        this->A = A;
        R_factor = std::make_shared<Matrix<m, n>>(*A);
        tau = std::make_shared<Matrix<m, 1>>(double(0.0));
        int res = LAPACKE_dgerqf(LAPACK_ROW_MAJOR, m, n, R_factor->getPtr(), n, tau->getPtr());
        if(res != 0)
        {
            LOGE << io::xprintf("Result is %d which is nonzero", res);
        }
    }

    template <uint m, uint n>
    std::shared_ptr<Matrix<m, m>> RQFactorization<m, n>::getRMatrix()
    {
        if(A != nullptr)
        {
            std::shared_ptr<Matrix<m, m>> R = std::make_shared<Matrix<m, m>>(0.0);
            for(uint32_t i = 0; i != m; i++) // Row
            {
                for(uint32_t j = i; j != m; j++) // Column filled only from i
                {
                    if((*R_factor)(j, j + n - m) < 0) // Constructing normal form in which I put
                                                      // nonnegative values to the diagonal
                    {
                        (*R)(i, j) = -(*R_factor)(i, j + n - m);
                    } else
                    {
                        (*R)(i, j) = (*R_factor)(i, j + n - m);
                    }
                }
            }
            return R;
        } else
        {
            io::throwerr(
                "You have to factorize matrix first by calling factor to get its decomposition.");
            return nullptr;
        }
    }

    template <uint m, uint n>
    std::shared_ptr<Matrix<m, n>> RQFactorization<m, n>::getQMatrix()
    {
        if(A != nullptr)
        {
            std::shared_ptr<Matrix<m, n>> Q = std::make_shared<Matrix<m, n>>(*R_factor);
            int res = LAPACKE_dorgrq(LAPACK_ROW_MAJOR, m, n, m, Q->getPtr(), n, tau->getConstPtr());
            if(res != 0)
            {
                LOGE << io::xprintf("Result is %d which is nonzero", res);
            }
            for(uint32_t i = 0; i != m; i++) // Row
            {
                if((*R_factor)(i, i + n - m) < 0) // Constructing normal form in which I put
                                                  // nonnegative values to the diagonal
                {
                    for(uint32_t j = 0; j != n; j++)
                    {
                        (*Q)(i, j) = -(*Q)(i, j);
                    }
                }
            }
            return Q;
        } else
        {
            io::throwerr(
                "You have to factorize matrix first by calling factor to get its decomposition.");
            return nullptr;
        }
    }

} // namespace matrix
} // namespace CTL
