#pragma once

namespace CTL::utils {
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

    double regularizeInverted(float s, double lambda)
    {
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
            double factor = (s * s) / (s * s + lambda * lambda);
            return factor / s;
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
        u = new float[n * n];
        s = new float[n];
        vt = new float[n * n];
        int info;
        info = LAPACKE_sgesdd(LAPACK_ROW_MAJOR, 'A', n, n, A, n, s, u, n, vt, n);
        // Everything will be stored row major, vt consecutive vectors, u vectors in columns non
        // consecutive
        if(info != 0)
        {
            io::throwerr("Info is %d which is nonzero value", info);
        }
        double lambda = lambda_rel * s[0];
        for(int i = 0; i != n * n; i++)
            A[i] = 0.0;
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
    }

private:
    double lambda_rel;
    bool truncate;
};

} // namespace CTL::utils
