#include "MATRIX/Matrix.hpp"

namespace CTL {
namespace matrix {

    Matrix::Matrix(uint32_t m, uint32_t n)
        : m(m)
        , n(n)
    {
        size = m * n;
        A = new double[size]();
    }

    Matrix::Matrix(uint32_t m, uint32_t n, const double x)
        : m(m)
        , n(n)
    {
        size = m * n;
        A = new double[size];
        std::fill_n(A, size, x);
    }

    Matrix::Matrix(uint32_t m, uint32_t n, const double* A_copy)
        : m(m)
        , n(n)
    {
        size = m * n;
        A = new double[size];
        std::copy(A_copy, A_copy + size, A);
    }

    Matrix::Matrix(uint32_t m, uint32_t n, std::initializer_list<double> A_copy)
        : m(m)
        , n(n)
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
    Matrix::~Matrix()
    {
        if(A != nullptr)
            delete[] A;
    }

    /// Copy constructor
    Matrix::Matrix(const Matrix& b)
        : Matrix(b.m, b.n, b.A)
    {
        // LOGD << "Caling Copy constructor of Matrix.";
    }

    /**Copy assignment
     *
     */
    Matrix& Matrix::operator=(const Matrix& b)
    {
        if(&b != this)
        {
            if(this->m == b.m && this->n == b.n)
            {
                memcpy(this->A, b.A, size * sizeof(double));
            } else
            {
                std::string msg = "Objects are incompatible and its dimensions differ. Can not use "
                                  "copy assignment operator!";
                LOGE << msg;
                throw new std::runtime_error(msg);
            }
        }
        return *this;

    } // copy assignment, tmp is to solve situation when assigning to itself

    Matrix::Matrix(Matrix&& b)
    {
        this->A = b.A;
        this->size = b.size;
        this->m = b.m;
        this->n = b.n;
        b.A = nullptr;
    } // Move constructor

    Matrix& Matrix::operator=(Matrix&& b)
    {
        if(&b != this)
        {
            delete[] this->A;
            this->A = b.A;
            this->size = b.size;
            this->m = b.m;
            this->n = b.n;
            b.A = nullptr;
        }
        return *this;
    } // Move assignment

    Matrix Matrix::minorSubMatrix(uint32_t i, uint32_t j) const
    {
        if(m == 0 || n == 0)
        {
            std::string errMsg
                = io::xprintf("Can not further subdivide matrix with dimension zero!", i, j);
            LOGE << errMsg;
            throw std::runtime_error(errMsg);
        }
        if(i > m - 1 || j > n - 1)
        {
            std::string errMsg = io::xprintf(
                "For minor submatrix you have to specify valid element (i, j) but not (%d, %d)!", i,
                j);
            LOGE << errMsg;
            throw std::runtime_error(errMsg);
        }
        Matrix out(m - 1, n - 1);
        for(uint32_t a = 0; a != m; a++)
        {
            for(uint32_t b = 0; b != n; b++)
            {
                if(i > a)
                {
                    if(j > b)
                    {
                        out(a, b) = (*this)(a, b);
                    } else if(j < b)
                    {
                        out(a, b - 1) = (*this)(a, b);
                    }
                } else if(i < a)
                {
                    if(j > b)
                    {
                        out(a - 1, b) = (*this)(a, b);
                    } else if(j < b)
                    {
                        out(a - 1, b - 1) = (*this)(a, b);
                    }
                }
            }
        }
        return (out);
    }

    Matrix Matrix::T() const
    {
        Matrix transposed(n, m);
        for(uint32_t i = 0; i < m; i++)
        {
            for(uint32_t j = 0; j < n; j++)
            {
                transposed(j, i) = this->get(i, j);
            }
        }
        return transposed;
    }

    /**Matrix multiplication is performed by BLAS
     *
     *Function dgemm is called.
     *https://software.intel.com/en-us/mkl-tutorial-c-multiplying-matrices-using-dgemm
     *https://software.intel.com/en-us/mkl-developer-reference-c-cblas-gemm
     */
    Matrix Matrix::operator*(const Matrix& B) const
    {
        Matrix out(this->m, B.n);
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, B.n, n, 1.0, getConstPtr(), n,
                    B.getConstPtr(), B.n, 0.0, out.getPtr(), B.n);
        return out;
    }

    Matrix Matrix::operator-(const Matrix& B) const
    {
        Matrix out(m, n);
        for(size_t i = 0; i != size; i++)
        {
            out.A[i] = A[i] - B.A[i];
        }
        return out;
    }

    std::string Matrix::toString(std::string name) const
    {
        std::ostringstream os;
        for(unsigned int i = 0; i != m; ++i)
        {
            if(i == 0)
            {
                os << io::xprintf("%s = |", name.c_str());
            } else
            {
                os << "    |";
            }
            for(unsigned int j = 0; j != n; ++j)
            {
                if(j != 0)
                    os << " ";
                os << std::setw(9) << std::fixed << std::setfill(' ') << std::setprecision(3)
                   << static_cast<double>((*this)(i, j));
            }
            os << "|";
            os << "\n";
        }
        os << std::endl;
        return os.str();
    }

    double Matrix::norm(uint16_t p) const
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
