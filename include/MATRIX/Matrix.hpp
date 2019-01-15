#pragma once

#include <plog/Log.h>

#include "mkl.h"

namespace CTL {
namespace matrix {
    /// Row major matrix with m rows and n columns
    template <uint m, uint n>
    class Matrix
    {
    public:
        /// Constructor creates matrix with given dimensions initialized to zero
        Matrix();
        /// Init matrix by single value
        Matrix(const double X);
        /// Constructor creates matrix with given dimensions initialized to zero
        Matrix(const double* A_copy);
        /// Constructor from initializer list
        Matrix(std::initializer_list<double> A_copy);
        /// Destructor
        ~Matrix();
        /// Copy constructor
        Matrix(const Matrix<m, n>& b);
        // Copy assignment
        Matrix<m, n>& operator=(const Matrix<m, n>& b);
        // Move constructor
        Matrix(Matrix<m, n>&& b);
        // Move assignment
        Matrix<m, n>& operator=(Matrix<m, n>&& other);

        static Matrix<m, n> unitDiagonal()
        {
            Matrix<m, n> X;
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

        /**ln norm of the matrix, p=0 means l\infty norm*/
        double norm(uint16_t p = 2) const;

        // Transpose matrix
        Matrix<n, m> T() const;

        // Matrix multiplication
        template <uint n2>
        Matrix<m, n2> operator*(const Matrix<n, n2>& B) const;
        Matrix<m, n> operator-(const Matrix<m, n>& rhs) const;

    private:
        double* A;
        size_t size;
    };

    template <uint m, uint n>
    Matrix<m, n>::Matrix()
    {
        size = m * n;
        A = new double[size]();
    }

    template <uint m, uint n>
    Matrix<m, n>::Matrix(const double* A_copy)
    {
        size = m * n;
        A = new double[size];
        std::copy(A_copy, A_copy + size, A);
    }

    template <uint m, uint n>
    Matrix<m, n>::Matrix(std::initializer_list<double> A_copy)
    {
        size = m * n;
        if(A_copy.size() != size)
        {
            io::throwerr("The class must be initialized with proper sized list of %d elements not "
                         "%d elements.",
                         size, A_copy.size());
        }
        A = new double[size];
        std::copy(A_copy.begin(), A_copy.end(), A);
    }

    // Guards stealed object removal
    template <uint m, uint n>
    Matrix<m, n>::~Matrix()
    {
        if(A != nullptr)
            delete[] A;
    }

    template <uint m, uint n>
    Matrix<m, n>::Matrix(const double x)
    {
        size = m * n;
        A = new double[size];
        std::fill_n(A, size, x);
    }

    /// Copy constructor
    template <uint m, uint n>
    Matrix<m, n>::Matrix(const Matrix<m, n>& b)
        : Matrix<m, n>(b.A)
    {
        LOGD << "Caling Copy constructor of Matrix.";
    }

    /**Copy assignment
     *
     */
    template <uint m, uint n>
    Matrix<m, n>& Matrix<m, n>::operator=(const Matrix<m, n>& b)
    {
        if(&b != this)
        {
            memcpy(this->A, b.A, size * sizeof(double));
        }
        return *this;

    } // copy assignment, tmp is to solve situation when assigning to itself

    template <uint m, uint n>
    Matrix<m, n>::Matrix(Matrix<m, n>&& b)
    {
        this->A = b.A;
        this->size = b.size;
        b.A = nullptr;
    } // Move constructor

    template <uint m, uint n>
    Matrix<m, n>& Matrix<m, n>::operator=(Matrix<m, n>&& b)
    {
        if(&b != this)
        {
            delete[] this->A;
            this->A = b.A;
            this->size = b.size;
            b.A = nullptr;
        }
        return *this;
    } // Move assignment

    template <uint m, uint n>
    Matrix<n, m> Matrix<m, n>::T() const
    {
        Matrix<n, m> ret;
        double* scrPtr = A;
        uint32_t row, column;

        for(row = 0; row < m; ++row)
            for(column = 0; column < n; ++column)
            {
                ret(column, row) = *scrPtr;
                ++scrPtr;
            }
        return ret;
    }

    /**Matrix multiplication is performed by BLAS
     *
     *Function dgemm is called.
     *https://software.intel.com/en-us/mkl-tutorial-c-multiplying-matrices-using-dgemm
     *https://software.intel.com/en-us/mkl-developer-reference-c-cblas-gemm
     */
    template <uint m, uint n>
    template <uint n2>
    Matrix<m, n2> Matrix<m, n>::operator*(const Matrix<n, n2>& B) const
    {
        Matrix<m, n2> out;
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n2, n, 1.0, getConstPtr(), n, B.getConstPtr(), n2,
                    0.0, out.getPtr(), n2);
        return out;
    }

    template <uint m, uint n>
    Matrix<m, n> Matrix<m, n>::operator-(const Matrix<m, n>& B) const
    {
        Matrix<m, n> out;
        for(size_t i = 0; i != size; i++)
        {
            out.A[i] = A[i] - B.A[i];
        }
        return out;
    }

    template <uint m, uint n>
    double Matrix<m, n>::norm(uint16_t p) const
    {
        if(p == 0)
        {
            double maximum = 0;
            for(size_t i = 0; i != size; i++)
            {
                maximum = std::max(maximum, std::abs(A[i]));
            }
            return maximum;
        } else
        {
            double ln = 0;
            double q = 1.0 / double(p);
            for(size_t i = 0; i != size; i++)
            {
                ln += std::pow(std::abs(A[i]), p);
            }
            return std::pow(ln, q);
        }
    }
} // namespace matrix
} // namespace CTL
