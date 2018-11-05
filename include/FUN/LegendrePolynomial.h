#pragma once

#include <cmath>

#include "stringFormatter.h"
#include "FUN/ScalarFunctionI.h"

namespace CTL {
namespace util {
    /**Interface that will be implemented to contain actual function supported on the interval
     * [start, end]
     *
     */
    class LegendrePolynomial : ScalarFunctionI
    {
    public:
        LegendrePolynomial(int degree, double start, double end)
            : ScalarFunctionI(start, end)
        {
            if(degree < 0)
            {
                io::throwerr("Polynomial degree must be greater then zero and not %d as supplied.",
                             degree);
            }
            this->degree = degree;
        } ///< Inits the function

        /**Provides value of the function at particular point.
         *
         *It is based on the formula
         *\f{equation}{
         *	P_n(x) = \sum_{k=0}^n \binom{n}{k} \binom{n+k}{k}(\frac{x-1}{2})^k
         *\f}
         */
        double valueAt(double t) const override
        {
            double x = transformToSupport(t);
            x = (x - 1) / 2;
            double val = 0.0;
            for(int k = 0; k != degree + 1; k++)
            {
                val += nChoosek(degree, k) * nChoosek(degree + k, k) * std::pow(x, k);
            }
            return val;
        }

        /**This one computes values of the polynomials starting from zeroth degree into the array.
         *
         */
        void valuesAt(double t, double* array)
        {
            // LOGW << io::xprintf("Degree is %d", degree);
            double x = transformToSupport(t);
            array[0] = 1.0;
            if(degree > 0)
            {
                array[1] = x;
            }
            if(degree > 1)
            {
                array[2] = ((3.0 * x * x) - 1.0) * 0.5;
            }
            for(int n = 3; n != degree + 1; n++)
            {
                array[n]
                    = (((2.0 * n) - 1.0) * x * array[n - 1] - (((double)n - 1.0) * array[n - 2]))
                    / (double)n;
            }
        }

    private:
        /**Function that transforms the value t on the interval [start, end] to the value t' on the
         * interval [-1,1] that is support of Legendre polynomials.
         *
         */
        double transformToSupport(double t) const { return (2 * (t - start) / (end - start)) - 1; }

        unsigned nChoosek(unsigned n, unsigned k) const
        {
            if(k > n)
                return 0;
            if(k * 2 > n)
                k = n - k;
            if(k == 0)
                return 1;

            int result = n;
            for(int i = 2; i <= k; ++i)
            {
                result *= (n - i + 1);
                result /= i;
            }
            return result;
        }

        inline double P0(double x) const { return 1.0; }

        // n = 1
        inline double P1(double x) const { return x; }

        // n = 2
        inline double P2(double x) const { return ((3.0 * x * x) - 1.0) * 0.5; }

        int degree;
    };

} // namespace util
} // namespace CTL
