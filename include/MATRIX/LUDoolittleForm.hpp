#pragma once
// Logging, testing
#include <plog/Log.h>

// External libraries
#include <cmath>
#include <memory>

// Internal libraries
#include "stringFormatter.h"
#include "MATRIX/Matrix.hpp"
#include "MATRIX/SquareMatrix.hpp"

namespace CTL {
namespace matrix {
    template <uint N>
    class LUDoolittleForm
    {
    public:
        LUDoolittleForm(std::shared_ptr<SquareMatrix<N>> LU,
                        std::shared_ptr<std::array<int, N>> P,
                        bool oddSwapParity)
            : LU(LU)
            , P(P)
            , oddSwapParity(oddSwapParity)
        {
        }

        double getDeterminant();
        double getLog10Determinant();
        SquareMatrix<N> inverseMatrix();
        Matrix<N, 1> solve(Matrix<N, 1> vector);
        SquareMatrix<N> getLMatrix();
        SquareMatrix<N> getUMatrix();
        SquareMatrix<N> getPAMatrix();
        std::array<int, N> getPermutation();
        bool getOddSwapParity();
        static LUDoolittleForm<N> LUDecomposeDoolittle(const SquareMatrix<N>& m,
                                                       double minimalPivotSize);

    private:
        std::shared_ptr<SquareMatrix<N>> LU;
        std::shared_ptr<std::array<int, N>> P;
        Matrix<N, 1> forwardSubstitute(Matrix<N, 1> b);
        Matrix<N, 1> backwardSubstitute(Matrix<N, 1> b);
        static int findPivot(const SquareMatrix<N>& M,
                             const std::array<int, N>& P,
                             int k,
                             double minimalPivotSize);
        bool oddSwapParity;
    };

    template <uint N>
    bool LUDoolittleForm<N>::getOddSwapParity()
    {
        return oddSwapParity;
    }

    template <uint N>
    SquareMatrix<N> LUDoolittleForm<N>::getLMatrix()
    {
        SquareMatrix<N> L;
        for(int i = 0; i != N; i++)
        {
            for(int j = 0; j != i; j++)
            {
                L(i, j) = (*LU)((*P)[i], j);
            }
            L(i, i) = 1.0;
            for(int j = i + 1; j != N; j++)
            {
                L(i, j) = 0.0;
            }
        }
        return L;
    }

    template <uint N>
    SquareMatrix<N> LUDoolittleForm<N>::getUMatrix()
    {
        SquareMatrix<N> U;
        for(int i = 0; i != N; i++)
        {
            for(int j = 0; j != i; j++)
            {
                U(i, j) = 0.0;
            }
            for(int j = i; j != N; j++)
            {
                U(i, j) = (*LU)((*P)[i], j);
            }
        }
        return U;
    }

    template <uint N>
    SquareMatrix<N> LUDoolittleForm<N>::getPAMatrix()
    {
	SquareMatrix<N> U = this->getUMatrix();
	SquareMatrix<N> L = this->getLMatrix();
        return L*U;
    }

    template <uint N>
    std::array<int, N> LUDoolittleForm<N>::getPermutation()
    {
        std::array<int, N> x;
        for(int i = 0; i != N; i++)
        {
            x[i] = (*P)[i];
        }
        return x;
    }

    /**Solution of the system A x = b, for b given.
     *
     * We are solving P A x = P b, L U x = P b, U x  = c = L^{-1} P b  by forward substitution, x =
     * U^{-1} c by back substitution
     */
    template <uint N>
    Matrix<N, 1> LUDoolittleForm<N>::solve(Matrix<N, 1> b)
    {
        // First we compute forward substitution
        Matrix<N, 1> c = forwardSubstitute(b);
        Matrix<N, 1> x = backwardSubstitute(c);
        return (x);
    }

    /**Forward substitution to solve LU decomposed linear system.
     *
     *See that the resulting vector is unpermutated, since its individual components multiply rows
     *of L in the equation Lc = Pb.
     *
     */
    template <uint N>
    Matrix<N, 1> LUDoolittleForm<N>::forwardSubstitute(Matrix<N, 1> b)
    {
        Matrix<N, 1> c;
        for(int i = 0; i != N; i++)
        {
            c(i, 0) = b((*P)[i], 0);
            for(int j = 0; j != i; j++)
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
    template <uint N>
    Matrix<N, 1> LUDoolittleForm<N>::backwardSubstitute(Matrix<N, 1> b)
    {
        Matrix<N, 1> x;
        for(int i = N - 1; i != -1; i--)
        {
            double divisor = (*LU)((*P)[i], i);
            x(i, 0) = b(i, 0);
            for(int j = i + 1; j != N; j++)
            {
                x(i, 0) -= (*LU)((*P)[i], j) * x(j, 0);
            }
            x(i, 0) /= divisor;
        }
        return x;
    }

    template <uint N>
    SquareMatrix<N> LUDoolittleForm<N>::inverseMatrix()
    {
        SquareMatrix<N> out;
        Matrix<N, 1> vec;
        for(int i = 0; i != N; i++)
        {
            Matrix<N, 1> eu(0.0);
            eu(i, 0) = 1;
            vec = solve(eu);
            //	LOGD << "Solution of " << eu.info() << "is" << vec.info();
            for(int j = 0; j != N; j++)
            {
                out(j, i) = vec(j,0);
            }
        }
        return out;
    }

    /**Compute determinant of the matrix based on LU decomposition.
     */
    template <uint N>
    double LUDoolittleForm<N>::getDeterminant()
    {
        double determinant = 1;
        for(int i = 0; i != N; i++)
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
    template <uint N>
    double LUDoolittleForm<N>::getLog10Determinant()
    {
        double determinant = 0;
        bool negativeSign = oddSwapParity;
        for(int i = 0; i != N; i++)
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
    template <uint N>
    LUDoolittleForm<N> LUDoolittleForm<N>::LUDecomposeDoolittle(const SquareMatrix<N>& M,
                                                                double minimalPivotSize)
    {
        int swapParity = 0;
        std::shared_ptr<SquareMatrix<N>> LU = std::make_shared<SquareMatrix<N>>(M);
        std::shared_ptr<std::array<int, N>> P = std::make_shared<std::array<int, N>>();
        for(int i = 0; i != N; i++)
        {
            (*P)[i]
                = i; // On i-th position is the row index in the original matrix that is on the i-th
                     // position in the PA matrix.
        }
        for(int k = 0; k != N; k++) // K-th step of the Gauss elimination
        {
            //    LOGD << io::string_format("Before %d-th step", k) << std::endl << LU->info();
            int p = LUDoolittleForm<N>::findPivot(*LU, *P, k, minimalPivotSize);
            if(p != k)
            {
                swapParity++;
                int tmp = (*P)[k];
                (*P)[k] = (*P)[p];
                (*P)[p] = tmp;
            }
            double pivot = (*LU)((*P)[k], k);
            //    LOGD << io::string_format("Pivot value is %f.", pivot);
            for(int i = k + 1; i != N; i++)
            {
                double l = (*LU)((*P)[i], k) / pivot;
                //       LOGD << io::string_format("For i=%d, j=%d, out[i,j]=%f, l=%f.",(*P)[i], k,
                //       (*LU)((*P)[i], k), l);
                (*LU)((*P)[i], k) = l;
                for(int j = k + 1; j != N; j++)
                {
                    (*LU)((*P)[i], j) -= (*LU)((*P)[k], j) * l;
                }
            }
            //    LOGD << io::string_format("After %d-th step", k) << std::endl << LU->info();
        }
        LUDoolittleForm<N> lu(LU, P, (swapParity % 2));
        return lu;
    }

    template <uint N>
    int LUDoolittleForm<N>::findPivot(const SquareMatrix<N>& M,
                                      const std::array<int, N>& P,
                                      int k,
                                      double minimalPivotSize)
    {
        double curMax = std::abs(M(P[k], k));
        int p = k;
        for(int i = k + 1; i != N; i++)
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
            LOGE << errMsg;
            throw std::runtime_error(errMsg);
        }
        return p;
    }

} // namespace util
} // namespace CTL
