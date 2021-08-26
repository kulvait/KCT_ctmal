#include "MATRIX/RQFactorization.hpp"

namespace KCT {
namespace matrix {
    void RQFactorization::factorize(std::shared_ptr<Matrix> A)
    {
        uint32_t m = A->dimm();
        uint32_t n = A->dimn();
        if(m > n)
        {
            io::throwerr("Initial matrix must be in the form m <= n but m=%d, n=%d.", m, n);
        }
        this->A = A;
        R_factor = std::make_shared<Matrix>(*A);
        tau = std::make_shared<Matrix>(m, 1, double(0.0));
        int res = LAPACKE_dgerqf(LAPACK_ROW_MAJOR, m, n, R_factor->getPtr(), n, tau->getPtr());
        if(res != 0)
        {
            LOGE << io::xprintf("Result is %d which is nonzero", res);
        }
    }

    std::shared_ptr<Matrix> RQFactorization::getRMatrix()
    {
        if(A != nullptr)
        {
            uint32_t m = A->dimm();
            uint32_t n = A->dimn();
            std::shared_ptr<Matrix> R = std::make_shared<Matrix>(m, m, 0.0);
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

    std::shared_ptr<Matrix> RQFactorization::getQMatrix()
    {
        if(A != nullptr)
        {
            uint32_t m = A->dimm();
            uint32_t n = A->dimn();
            std::shared_ptr<Matrix> Q = std::make_shared<Matrix>(*R_factor);
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
} // namespace KCT
