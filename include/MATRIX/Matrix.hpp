#pragma once

#include <plog/Log.h>

// STL
#include <cmath>

// External libraries
#include "mkl.h"
#include "stringFormatter.h"

namespace KCT {
namespace matrix {
    /// Row major matrix with m rows and n columns
    class Matrix
    {
    public:
        /// Constructor creates matrix with given dimensions initialized to zero
        Matrix(uint32_t m, uint32_t n);
        /// Init matrix by single value
        Matrix(uint32_t m, uint32_t n, const double X);
        /// Constructor creates matrix with given dimensions initialized to zero
        Matrix(uint32_t m, uint32_t n, const double* A_copy);
        /// Constructor from initializer list
        Matrix(uint32_t m, uint32_t n, std::initializer_list<double> A_copy);
        /// Destructor
        ~Matrix();
        /// Copy constructor
        Matrix(const Matrix& b);
        // Copy assignment
        Matrix& operator=(const Matrix& b);
        // Move constructor
        Matrix(Matrix&& b);
        // Move assignment
        Matrix& operator=(Matrix&& other);

/**
* @brief Create submatrix formed by removal of the row and column coresponding to m and n.
*
* @param m
* @param n
*
* @return 
*/
	Matrix minorSubMatrix(uint32_t m, uint32_t n) const;

        static Matrix unitDiagonal(uint32_t m, uint32_t n)
        {
            Matrix X(m, n);
            for(uint32_t i = 0; i != std::min(m, n); i++)
            {
                X(i, i) = 1.0;
            }
            return X;
        }

        /**Get data pointer.
         *
         */
        double* getPtr() { return A; }

        /**Get const data pointer.
         *
         */
        const double* getConstPtr() const { return A; }

        /**Unchecked element acess.
         *
         */
        double& operator()(uint32_t i, uint32_t j) { return A[i * n + j]; }
        double operator()(uint32_t i, uint32_t j) const { return A[i * n + j]; }
        double get(uint32_t i, uint32_t j) const { return A[i * n + j]; };

        /**ln norm of the matrix, p=0 means l\infty norm*/
        double norm(uint16_t p = 2) const;
        uint32_t dimm() const { return m; };
        uint32_t dimn() const { return n; };

        // Transpose matrix
        Matrix T() const;

        /**Creates string representation of the matrix
         *
         *Do not use for large matrices.
         */
        std::string toString(std::string name = "A") const;
        // Matrix multiplication
        Matrix operator*(const Matrix& B) const;
        Matrix operator-(const Matrix& rhs) const;
        friend Matrix operator*(const double& alpha, const Matrix& rhs)
        {
            Matrix out(rhs.m, rhs.n);
            for(size_t i = 0; i != rhs.size; i++)
            {
                out.A[i] = alpha * rhs.A[i];
            }
            return out;
        }

    protected:
        double* A;
        size_t size;
        uint32_t m, n;
    };

} // namespace matrix
} // namespace KCT
