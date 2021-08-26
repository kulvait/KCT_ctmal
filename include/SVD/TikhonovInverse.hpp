#pragma once

#include "mkl_lapacke.h" //Lapacke

namespace KCT::utils {
/**Construct approximate inverse of square matrix using SVD, where Tikhonov regularization is
 * utilized.
 *
 */
class TikhonovInverse
{
public:
    /** When truncate is true computes truncated SVD.
     */
    TikhonovInverse(double lambda_rel, bool truncate)
    {
        this->lambda_rel = lambda_rel;
        this->truncate = truncate;
    }

    /** Precompute convolution matrix in the format as on page 10, Fieselmann at. al, 2011.
     *
     *@param[in] granularity Number of elements in aif sequence, A have to be prealocated to have
     *granularity*granularity elements.
     *@param[in] aif Values of arthery input function to make deconvolution matrix from.
     *@param[out] A The matrix prealocated to size granularity*granularity to be filled by
     *deconvolution data, row major alignment.
     *
     */
    static void precomputeConvolutionMatrix(const uint32_t granularity, const float* aif, float* A)
    {
        for(uint32_t i = 0; i != granularity; i++)
            for(uint32_t j = 0; j != granularity; j++)
            {
                if(j > i)
                {
                    A[i * granularity + j] = 0.0;
                } else
                {
                    A[i * granularity + j] = aif[i - j];
                }
            }
    }

    double regularizeInverted(float s, double lambda)
    {
        if(s == 0)
        {
            return 0.0;
        }
        if(truncate)
        {
            if(s < lambda)
            {
                return 0;
            } else
            {
                return 1.0 / s;
            }
        } else
        {
            double factor = s / (s * s + lambda * lambda);
            return factor;
        }
    }

    /**
     *
     * A is the matrix in row major order
     * n is number of rows of the matrix such that matrix has n*n elements allocated
     * in array A is outputted the result of the computation
     */
    void computePseudoinverse(float* A, int n)
    {
        float *u, *s, *vt;
        float* superb = new float[n - 1];
        u = new float[n * n];
        s = new float[n];
        vt = new float[n * n];
        int info;
        info = LAPACKE_sgesvd(LAPACK_ROW_MAJOR, 'A', 'A', n, n, A, n, s, u, n, vt, n, superb);
        // Everything will be stored row major, vt consecutive vectors, u vectors in columns non
        // consecutive
        if(info != 0)
        {
            io::throwerr("Info is %d which is nonzero value", info);
        }
        double lambda = lambda_rel * s[0];
        std::fill(A, A + n * n, float(0.0));
        for(int k = 0; k != n; k++)
        {
            double sigma = regularizeInverted(s[k], lambda);
            if(sigma != 0.0)
            {
                for(int i = 0; i != n; i++)
                {
                    for(int j = 0; j != n; j++)
                    {
                        A[i * n + j] += vt[k * n + i] * sigma * u[n * j + k];
                    }
                }
            }
        }
        delete[] u;
        delete[] s;
        delete[] vt;
        delete[] superb;
    }

private:
    double lambda_rel;
    bool truncate;
};

} // namespace KCT::utils
