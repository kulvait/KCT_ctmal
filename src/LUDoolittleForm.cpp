#include "MATRIX/LUDoolittleForm.hpp"

namespace CTL {
namespace matrix {
    bool LUDoolittleForm::getOddSwapParity() { return oddSwapParity; }

    SquareMatrix LUDoolittleForm::getLMatrix()
    {
        SquareMatrix L(n);
        for(uint32_t i = 0; i != n; i++)
        {
            for(uint32_t j = 0; j != i; j++)
            {
                L(i, j) = (*LU)((*P)[i], j);
            }
            L(i, i) = 1.0;
            for(uint32_t j = i + 1; j != n; j++)
            {
                L(i, j) = 0.0;
            }
        }
        return L;
    }

    SquareMatrix LUDoolittleForm::getUMatrix()
    {
        SquareMatrix U(n);
        for(uint32_t i = 0; i != n; i++)
        {
            for(uint32_t j = 0; j != i; j++)
            {
                U(i, j) = 0.0;
            }
            for(uint32_t j = i; j != n; j++)
            {
                U(i, j) = (*LU)((*P)[i], j);
            }
        }
        return U;
    }

    SquareMatrix LUDoolittleForm::getPAMatrix()
    {
        SquareMatrix U = this->getUMatrix();
        SquareMatrix L = this->getLMatrix();
        return L * U;
    }

    std::vector<uint32_t> LUDoolittleForm::getPermutation()
    {
        std::vector<uint32_t> x = *P;
        return x; // Deep copy
    }

    /**Solution of the system A x = b, for b given.
     *
     * We are solving P A x = P b, L U x = P b, U x  = c = L^{-1} P b  by forward substitution, x =
     * U^{-1} c by back substitution
     */
    Matrix LUDoolittleForm::solve(Matrix b)
    {
        // First we compute forward substitution
        Matrix c = forwardSubstitute(b);
        Matrix x = backwardSubstitute(c);
        return (x);
    }

    /**Forward substitution to solve LU decomposed linear system.
     *
     *See that the resulting vector is unpermutated, since its individual components multiply rows
     *of L in the equation Lc = Pb.
     *
     */
    Matrix LUDoolittleForm::forwardSubstitute(Matrix b)
    {
        Matrix c(n, 1);
        for(uint32_t i = 0; i != n; i++)
        {
            c(i, 0) = b((*P)[i], 0);
            for(uint32_t j = 0; j != i; j++)
            {
                c(i, 0) -= (*LU)((*P)[i], j) * c(j, 0);
            }
        }
        return c;
    }

    /**Backward substitution to solve LU decomposed linear system.
     *
     *See that the resulting vector is unpermutated since its individual components multiply rows of
     *U x = b, in this case b is also unpermutated.
     *
     */
    Matrix LUDoolittleForm::backwardSubstitute(Matrix b)
    {
        Matrix x(n, 1);
        for(int i = n - 1; i != -1; i--)
        {
            double divisor = (*LU)((*P)[i], i);
            x(i, 0) = b(i, 0);
            for(int j = i + 1; j < int(n); j++)
            {
                x(i, 0) -= (*LU)((*P)[i], j) * x(j, 0);
            }
            x(i, 0) /= divisor;
        }
        return x;
    }

    SquareMatrix LUDoolittleForm::inverseMatrix()
    {
        SquareMatrix out(n);
        Matrix vec(n, 1);
        for(uint32_t i = 0; i != n; i++)
        {
            Matrix eu(n, 1, 0.0);
            eu(i, 0) = 1;
            vec = solve(eu);
            //	LOGD << "Solution of " << eu.info() << "is" << vec.info();
            for(uint32_t j = 0; j != n; j++)
            {
                out(j, i) = vec(j, 0);
            }
        }
        return out;
    }

    /**Compute determinant of the matrix based on LU decomposition.
     */
    double LUDoolittleForm::getDeterminant()
    {
        double determinant = 1;
        for(uint32_t i = 0; i != n; i++)
        {
            determinant *= (*LU)((*P)[i], i);
        }
        if(oddSwapParity)
        {
            return -determinant;
        } else
        {
            return determinant;
        }
    }

    /* If the determinant is positive, it computes its log10 if its negative it computes -log10 of
     * the absolute value of the determinant.
     */
    double LUDoolittleForm::getLog10Determinant()
    {
        double determinant = 0;
        bool negativeSign = oddSwapParity;
        for(uint32_t i = 0; i != n; i++)
        {
            double val = (*LU)((*P)[i], i);
            if(val < 0)
            {
                negativeSign = !negativeSign;
            }
            determinant += std::log10(std::abs(val));
        }
        if(negativeSign)
        {
            return -determinant;
        } else
        {
            return determinant;
        }
    }

    /**LUP decomposition of a matrix, Doolittle factorization
     *
     *minimalPivotSize ... the minimal absolute value of the pivot to proceed the method
     *If there is a situation when there is no pivot of sufficient size, that is minimalPivotSize,
     *exception is thrown. Pivot in k-th step is choosen from the elements A^k_{ik}, i>=k of the
     *matrix A^k that is temporary matrix of the Gauss elimination. This methods performs Doolittle
     *factorization with partial pivoting, that means L is a lower triangular matrix with all 1 on
     *the diagonal. U is a upper triangular matrix. We transform matrix to the form PA = LU. P is a
     *permutation of matrix rows.
     */
    LUDoolittleForm LUDoolittleForm::LUDecomposeDoolittle(const SquareMatrix& M,
                                                          double minimalPivotSize)
    {
        int swapParity = 0;
        uint32_t n = M.dimn();
        std::shared_ptr<SquareMatrix> LU = std::make_shared<SquareMatrix>(M);
        std::shared_ptr<std::vector<uint32_t>> P = std::make_shared<std::vector<uint32_t>>();
        for(uint32_t i = 0; i != n; i++)
        {
            P->push_back(i);
            // On i-th position is the row index in the original matrix that is on the i-th
            // position in the PA matrix.
        }
        for(uint32_t k = 0; k != n; k++) // K-th step of the Gauss elimination
        {
            //    LOGD << io::string_format("Before %d-th step", k) << std::endl << LU->info();
            uint32_t p = LUDoolittleForm::findPivot(*LU, *P, k, minimalPivotSize);
            if(p != k)
            {
                swapParity++;
                int tmp = (*P)[k];
                (*P)[k] = (*P)[p];
                (*P)[p] = tmp;
            }
            double pivot = (*LU)((*P)[k], k);
            //    LOGD << io::string_format("Pivot value is %f.", pivot);
            for(uint32_t i = k + 1; i != n; i++)
            {
                double l = (*LU)((*P)[i], k) / pivot;
                //       LOGD << io::string_format("For i=%d, j=%d, out[i,j]=%f, l=%f.",(*P)[i], k,
                //       (*LU)((*P)[i], k), l);
                (*LU)((*P)[i], k) = l;
                for(uint32_t j = k + 1; j != n; j++)
                {
                    (*LU)((*P)[i], j) -= (*LU)((*P)[k], j) * l;
                }
            }
            //    LOGD << io::string_format("After %d-th step", k) << std::endl << LU->info();
        }
        LUDoolittleForm lu(LU, P, (swapParity % 2));
        return lu;
    }

    uint32_t LUDoolittleForm::findPivot(const SquareMatrix& M,
                                        const std::vector<uint32_t>& P,
                                        uint32_t k,
                                        double minimalPivotSize)
    {
        int n = M.dimn();
        double curMax = std::abs(M(P[k], k));
        int p = k;
        for(int i = k + 1; i != n; i++)
        {
            if(curMax < std::abs(M(P[i], k)))
            {
                p = i;
                curMax = std::abs(M(P[i], k));
            }
        }
        if(curMax <= minimalPivotSize)
        {
            std::string errMsg = io::xprintf(
                "Matrix is so close to singular that I can not find pivot of sufficient size.");
            //            LOGE << errMsg; ... here is not necessery to print msg when I use it for
            //            determining determinant even for singular matrices
            throw PivotingException(errMsg);
        }
        return p;
    }

    double LUDoolittleForm::determinant(const SquareMatrix& A)
    {
        if(A.dimm() == 0)
        {
            std::string msg = "Can not compute determinant of matrix with dimensions 0x0!";
            LOGE << msg;
            throw std::runtime_error(msg);
        }
        if(A.dimm() == 1)
        {
            return A.get(0, 0);
        }
        try
        {
            LUDoolittleForm lu = LUDoolittleForm::LUDecomposeDoolittle(A, 0.000001);
            return lu.getDeterminant();
        } catch(const matrix::PivotingException e)
        { // Pivoting failed so that determinant is close to zero, attempt to compute it anyway
            double det = 0.0;
            for(uint32_t j = 0; j != A.dimn(); j++)
            {
                Matrix minorSub = A.minorSubMatrix(0, j);
                SquareMatrix minorSubMatrix(minorSub);
                double minor = determinant(minorSubMatrix);
                if(j % 2 == 0)
                {
                    det += A(0, j) * minor;
                } else
                {
                    det -= A(0, j) * minor;
                }
            }
            return det;
        }
    }

} // namespace matrix
} // namespace CTL
