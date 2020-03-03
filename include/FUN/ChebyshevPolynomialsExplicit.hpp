#pragma once

#include <cmath>
#include <iostream>
#include <mutex>

#include "VectorFunctionI.h"
#include "stringFormatter.h"

namespace CTL {
namespace util {
    /** Class for evaluation Chebyshev polynomials of the first kind stretched from [-1,1] to a
     domain [start, end].
     *
     *Values at particular point of the domain are evaluated for Chebyshev polynomials. Values are
     *written to the array. Polynomials start from the polynomialDegree startReportingDegree and
     *ends by the polynomial of the polynomialDegree. polynomialDegree ... polynomialDegree of
     *polynomial to report last in the resulting vector
     *
     * For the definition see T_i(x) definition from
     *     https://en.wikipedia.org/wiki/Chebyshev_polynomials
     */
    class ChebyshevPolynomialsExplicit : public VectorFunctionI
    {
    public:
        ChebyshevPolynomialsExplicit(uint32_t polynomialDegree,
                                     double start,
                                     double end,
                                     bool constantOutsideInterval = true,
                                     uint32_t startReportingDegree = 0)
            : VectorFunctionI(polynomialDegree + 1 - startReportingDegree, start, end)
            , transformationSlope(double(2.0 / (end - start)))
            , transformationIntercept(-double((end + start) / (end - start)))
            , polynomialDegree(polynomialDegree)
            , constantOutsideInterval(constantOutsideInterval)
            , startReportingDegree(startReportingDegree)
        {
            if(startReportingDegree > polynomialDegree)
            {
                std::string msg = io::xprintf(
                    "Parameter sartReportingDegree=%d but it should be less or equal than "
                    "polynomialDegree=%d, correct!",
                    startReportingDegree, polynomialDegree);
                LOGD << msg;
                throw std::runtime_error(msg);
            }
            // Now precompute the values of chebyshev polynomials
            chebyshevCoefficientsD = new double[(polynomialDegree + 1) * (polynomialDegree + 1)];
            chebyshevCoefficientsF = new float[(polynomialDegree + 1) * (polynomialDegree + 1)];
            xnF = new float[polynomialDegree + 1];
            xnD = new double[polynomialDegree + 1];
            fillChebyshevCoefficients(chebyshevCoefficientsD, polynomialDegree);
            fillChebyshevCoefficients(chebyshevCoefficientsF, polynomialDegree);
        } ///< Inits the function

        ~ChebyshevPolynomialsExplicit()
        {
            if(chebyshevCoefficientsD != nullptr)
            {
                delete[] chebyshevCoefficientsD;
            }
            if(chebyshevCoefficientsF != nullptr)
            {
                delete[] chebyshevCoefficientsF;
            }
            if(xnF != nullptr)
            {
                delete[] xnF;
            }
            if(xnD != nullptr)
            {
                delete[] xnD;
            }
            chebyshevCoefficientsD = nullptr;
            chebyshevCoefficientsF = nullptr;
            xnF = nullptr;
            xnD = nullptr;
        }

        ChebyshevPolynomialsExplicit(const ChebyshevPolynomialsExplicit& other)
            : ChebyshevPolynomialsExplicit(
                  other.polynomialDegree, other.start, other.end, other.startReportingDegree)
        {

        } // Copy constructor

        ChebyshevPolynomialsExplicit& operator=(ChebyshevPolynomialsExplicit other)
        {
            if(this != &other)
            {
                VectorFunctionI::operator=(other);
                double interval = other.end - other.start;
                this->transformationSlope = double(2.0 / interval);
                this->transformationIntercept = -double((other.end + other.start) / interval);
                this->startReportingDegree = other.startReportingDegree;
                if(this->polynomialDegree != other.polynomialDegree)
                {
                    this->polynomialDegree = other.polynomialDegree;
                    if(chebyshevCoefficientsD != nullptr)
                    {
                        delete[] chebyshevCoefficientsD;
                    }
                    if(chebyshevCoefficientsF != nullptr)
                    {
                        delete[] chebyshevCoefficientsF;
                    }
                    if(xnF != nullptr)
                    {
                        delete[] xnF;
                    }
                    if(xnD != nullptr)
                    {
                        delete[] xnD;
                    }
                    chebyshevCoefficientsD = nullptr;
                    chebyshevCoefficientsF = nullptr;
                    xnF = nullptr;
                    xnD = nullptr;
                    xnF = new float[polynomialDegree + 1];
                    xnD = new double[polynomialDegree + 1];
                    chebyshevCoefficientsD
                        = new double[(this->polynomialDegree + 1) * (this->polynomialDegree + 1)];
                    chebyshevCoefficientsF
                        = new float[(this->polynomialDegree + 1) * (this->polynomialDegree + 1)];
                    fillChebyshevCoefficients(chebyshevCoefficientsD, polynomialDegree);
                    fillChebyshevCoefficients(chebyshevCoefficientsF, polynomialDegree);
                }
            }
            return *this;
        } // Assignment

#if DEBUG
        void plotFunctions(uint32_t granularity = 100,
                           std::shared_ptr<std::vector<std::string>> names = nullptr) override
        {
            if(names == nullptr)
            {
                names = std::make_shared<std::vector<std::string>>();
                for(uint32_t i = startReportingDegree; i <= polynomialDegree; ++i)
                {
                    names->push_back(io::xprintf("Chebyshev %d", i));
                }
            }
            VectorFunctionI::plotFunctions(granularity, names);
        }
#endif
        /*Construct the Chebyshev polynomials using polynomial basis.
         *
         * It is based on the formulas
         *
         *\f{equation}{
         *       P_0 (x) = 1, \quad P_1(x) = x,
         *\f}
         *
         *\f{equation}{
         *  P_n = 2 x P_{n-1}(x)-P_{n-2}.
         *\f}
         */
        static void fillChebyshevCoefficients(double* c, int deg)
        {
            std::fill(c, &c[(deg + 1) * (deg + 1)], double(0.0));
            c[0] = double(1.0);
            if(deg > 0)
            {
                c[deg + 1 + 1] = double(1.0);
            }
            if(deg > 1)
            {
                double factor_1 = 2.0, factor_2 = -1.0;
                for(int n = 2; n != deg + 1; n++)
                {
                    for(int i = 0; i != n; i++)
                    {
                        c[n * (deg + 1) + 1 + i] = factor_1 * c[(n - 1) * (deg + 1) + i];
                        c[n * (deg + 1) + i] += factor_2 * c[(n - 2) * (deg + 1) + i];
                    }
                }
            }
        }

        /*Construct the Chebyshev polynomials using polynomial basis.
         *
         * It is based on the formulas
         *
         *\f{equation}{
         *       P_0 (x) = 1, \quad P_1(x) = x,
         *\f}
         *
         *\f{equation}{
         *  P_n = \frac{(2n-1) x P_{n-1}(x)-(n-1)P_{n-2}}{n}.
         *\f}
         */
        static void fillChebyshevCoefficients(float* c, int deg)
        {
            std::fill(c, &c[(deg + 1) * (deg + 1)], float(0.0));
            c[0] = float(1.0);
            if(deg > 0)
            {
                c[deg + 1 + 1] = float(1.0);
            }
            if(deg > 1)
            {
                double factor_1 = 2.0f, factor_2 = -1.0f;
                for(int n = 2; n != deg + 1; n++)
                {
                    for(int i = 0; i != n; i++)
                    {
                        c[n * (deg + 1) + 1 + i] = factor_1 * c[(n - 1) * (deg + 1) + i];
                        c[n * (deg + 1) + i] += factor_2 * c[(n - 2) * (deg + 1) + i];
                    }
                }
            }
        }

        /**Values of the Chebyshev polynomials at specific time point.
         *
         * This function needs to be protected by mutex due to filling array of powers.
         * Alternative implementation is to create local array of powers.
         */
        void valuesAt(double t, double* array) const override
        {
            if(t < start)
            {
                if(constantOutsideInterval)
                {
                    valuesOutsideInterval(array);
                    return;
                }
                t = start;
            }
            if(t > end)
            {
                if(constantOutsideInterval)
                {
                    valuesOutsideInterval(array);
                    return;
                }
                t = end;
            }
            // std::lock_guard<std::mutex> guard(powerProtectionMutex);//Big overhead
            double x = transformToSupport(t);
            double* xn = new double[polynomialDegree + 1];
            xn[0] = 1.0;
            for(uint32_t n = 1; n < polynomialDegree + 1; n++)
            {
                // xnD[n] = std::pow(x, double(n)); // Numerically more stable but slower
                xn[n] = xn[n - 1] * x;
            }
            std::fill(array, &array[polynomialDegree - startReportingDegree + 1], double(0.0));
            for(uint32_t i = startReportingDegree; i < polynomialDegree + 1; i++)
            {
                for(uint32_t n = 0; n != i + 1; n++)
                {
                    array[i - startReportingDegree]
                        += chebyshevCoefficientsD[i * (polynomialDegree + 1) + n] * xn[n];
                }
            }
            delete[] xn;
        }

        /**Values of the Chebyshev polynomials at specific time point.
         *
         */
        void valuesAt(double t, float* array) const override
        {
            if(t < start)
            {
                if(constantOutsideInterval)
                {
                    valuesOutsideInterval(array);
                    return;
                }
                t = start;
            }
            if(t > end)
            {
                if(constantOutsideInterval)
                {
                    valuesOutsideInterval(array);
                    return;
                }
                t = end;
            }
            // std::lock_guard<std::mutex> guard(powerProtectionMutex);//Big overhead
            double x = transformToSupport(t);
            double xk = 1.0;
            float* xn = new float[polynomialDegree + 1];
            for(uint32_t n = 0; n < polynomialDegree + 1; n++)
            {
                // xnF[n] = std::pow(x, float(n)); // Numerically more stable but slower
                xn[n] = float(
                    xk); // Power is double precission and particular value is single precision
                xk *= x;
            }
            std::fill(array, &array[polynomialDegree - startReportingDegree + 1], float(0.0));
            for(uint32_t i = startReportingDegree; i < polynomialDegree + 1; i++)
            {
                for(uint32_t n = 0; n != i + 1; n++)
                {
                    array[i - startReportingDegree]
                        += chebyshevCoefficientsF[i * (polynomialDegree + 1) + n] * xn[n];
                }
            }
            delete[] xn;
        }

        /**Into the array insert the polynomial basis of Chebyshev polynomial of certain
         * polynomialDegree.
         *
         * array must be prealocated to the size deg+1
         * deg must be between 0 and polynomialDegree including
         */
        void polynomialValues(uint32_t deg, float* array)
        {
            if(deg > polynomialDegree)
            {
                io::throwerr("Degree %d must be in [0, %d]!", deg, polynomialDegree);
            }
            std::copy(&chebyshevCoefficientsF[deg * (polynomialDegree + 1)],
                      &chebyshevCoefficientsF[deg * (polynomialDegree + 1) + deg + 1], array);
        }

        /**Into the array insert the polynomial basis of Chebyshev polynomial of certain
         * polynomialDegree.
         *
         * array must be prealocated to the size deg+1
         * deg must be between 0 and polynomialDegree including
         */
        void polynomialValues(uint32_t deg, double* array)
        {
            if(deg > polynomialDegree)
            {
                io::throwerr("Degree %d must be in [0, %d]!", deg, polynomialDegree);
            }
            std::copy(&chebyshevCoefficientsD[deg * (polynomialDegree + 1)],
                      &chebyshevCoefficientsD[deg * (polynomialDegree + 1) + deg + 1], array);
        }

        void printChebyshevPolynomial(uint32_t deg)
        {
            if(deg > polynomialDegree)
            {
                io::throwerr("Degree %d must be in [0, %d]!", deg, polynomialDegree);
            }
            std::cout << io::xprintf("T_%d(x) = ", deg);
            if(chebyshevCoefficientsD[deg * (polynomialDegree + 1)] != 0)
            {
                std::cout << io::xprintf(" %.1f",
                                         chebyshevCoefficientsD[deg * (polynomialDegree + 1)]);
            }
            for(uint32_t i = 1; i < deg + 1; i++)
            {
                if(chebyshevCoefficientsD[deg * (polynomialDegree + 1) + i] != 0)
                {
                    double c = chebyshevCoefficientsD[deg * (polynomialDegree + 1) + i];
                    if(c < 0)
                    {
                        std::cout << io::xprintf(
                            " %.1fx^%d", chebyshevCoefficientsD[deg * (polynomialDegree + 1) + i],
                            i);
                    } else
                    {
                        std::cout << io::xprintf(
                            " +%.1fx^%d", chebyshevCoefficientsD[deg * (polynomialDegree + 1) + i],
                            i);
                    }
                }
            }
            std::cout << std::endl;
        }

        // Transformation to the interval [start,end] => [-1,1]
        double transformToSupport(double t) const
        {
            return (transformationSlope * t) + transformationIntercept;
        }

    private:
        void valuesOutsideInterval(float* array) const
        {
            for(uint32_t i = startReportingDegree; i < polynomialDegree + 1; i++)
            {
                if(i == 0)
                {
                    array[i - startReportingDegree] = 1.0f;
                } else
                {
                    array[i - startReportingDegree] = 0.0f;
                }
            }
        }

        void valuesOutsideInterval(double* array) const
        {
            for(uint32_t i = startReportingDegree; i < polynomialDegree + 1; i++)
            {
                if(i == 0)
                {
                    array[i - startReportingDegree] = 1.0;
                } else
                {
                    array[i - startReportingDegree] = 0.0;
                }
            }
        }
        /**Function that transforms the value t on the interval [start, end] to the value t' on the
         * interval [-1,1] that is support of Chebyshev polynomials.
         *
         */
        double transformationSlope;
        double transformationIntercept;
        uint32_t polynomialDegree;
        bool constantOutsideInterval;
        uint32_t startReportingDegree;

        double* chebyshevCoefficientsD;
        float* chebyshevCoefficientsF;
        double* xnD;
        float* xnF;
        mutable std::mutex powerProtectionMutex;
    };

} // namespace util
} // namespace CTL
