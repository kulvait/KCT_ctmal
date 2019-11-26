#pragma once

#include <cmath>
#include <iostream>
#include <mutex>

#include "VectorFunctionI.h"
#include "stringFormatter.h"

namespace CTL {
namespace util {
    /** Class for evaluation Legendre polynomials stretched from [-1,1] to a domain [start, end].
     *
     *Values at particular point of the domain are evaluated for Legendre polynomials. Values are
     *written to the array. Polynomials start from the polynomialDegree startReportingDegree and
     *ends by the polynomial of the polynomialDegree. polynomialDegree ... polynomialDegree of
     *polynomial to report last in the resulting vector
     */
    class LegendrePolynomialsExplicit : public VectorFunctionI
    {
    public:
        LegendrePolynomialsExplicit(uint32_t polynomialDegree,
                                    double start,
                                    double end,
                                    uint32_t startReportingDegree=0)
            : VectorFunctionI(polynomialDegree + 1 - startReportingDegree, start, end)
            , transformationSlope(double(2.0 / (end - start)))
            , transformationIntercept(-double((end + start) / (end - start)))
            , polynomialDegree(polynomialDegree)
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
            // Now precompute the values of legendre polynomials
            legendreCoefficientsD = new double[(polynomialDegree + 1) * (polynomialDegree + 1)];
            legendreCoefficientsF = new float[(polynomialDegree + 1) * (polynomialDegree + 1)];
            xnF = new float[polynomialDegree + 1];
            xnD = new double[polynomialDegree + 1];
            fillLegendreCoefficients(legendreCoefficientsD, polynomialDegree);
            fillLegendreCoefficients(legendreCoefficientsF, polynomialDegree);
        } ///< Inits the function

        ~LegendrePolynomialsExplicit()
        {
            if(legendreCoefficientsD != nullptr)
            {
                delete[] legendreCoefficientsD;
            }
            if(legendreCoefficientsF != nullptr)
            {
                delete[] legendreCoefficientsF;
            }
            if(xnF != nullptr)
            {
                delete[] xnF;
            }
            if(xnD != nullptr)
            {
                delete[] xnD;
            }
            legendreCoefficientsD = nullptr;
            legendreCoefficientsF = nullptr;
            xnF = nullptr;
            xnD = nullptr;
        }

        LegendrePolynomialsExplicit(const LegendrePolynomialsExplicit& other)
            : LegendrePolynomialsExplicit(
                  other.polynomialDegree, other.start, other.end, other.startReportingDegree)
        {

        } // Copy constructor

        LegendrePolynomialsExplicit& operator=(LegendrePolynomialsExplicit other)
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
                    if(legendreCoefficientsD != nullptr)
                    {
                        delete[] legendreCoefficientsD;
                    }
                    if(legendreCoefficientsF != nullptr)
                    {
                        delete[] legendreCoefficientsF;
                    }
                    if(xnF != nullptr)
                    {
                        delete[] xnF;
                    }
                    if(xnD != nullptr)
                    {
                        delete[] xnD;
                    }
                    legendreCoefficientsD = nullptr;
                    legendreCoefficientsF = nullptr;
                    xnF = nullptr;
                    xnD = nullptr;
                    xnF = new float[polynomialDegree + 1];
                    xnD = new double[polynomialDegree + 1];
                    legendreCoefficientsD
                        = new double[(this->polynomialDegree + 1) * (this->polynomialDegree + 1)];
                    legendreCoefficientsF
                        = new float[(this->polynomialDegree + 1) * (this->polynomialDegree + 1)];
                    fillLegendreCoefficients(legendreCoefficientsD, polynomialDegree);
                    fillLegendreCoefficients(legendreCoefficientsF, polynomialDegree);
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
                for(int i = startReportingDegree; i <= polynomialDegree; ++i)
                {
                    names->push_back(io::xprintf("Legendre %d", i));
                }
            }
            VectorFunctionI::plotFunctions(granularity, names);
        }
#endif
        /*Construct the Legendre polynomials using polynomial basis.
         *
         * It is based on the formulas
         *
         *\f{equation}{
         *       P_0 (x) = 1, \quad P_1(x) = x,
         *\f}
         *
         *\f{equation}{
         *  P_n = \frac{(2n-1) P_{n-1}(x)-(n-1)P_{n-2}}{n}.
         *\f}
         */
        static void fillLegendreCoefficients(double* c, int deg)
        {
            std::fill(c, &c[(deg + 1) * (deg + 1)], double(0.0));
            c[0] = double(1.0);
            if(deg > 0)
            {
                c[deg + 1 + 1] = double(1.0);
            }
            if(deg > 1)
            {
                c[2 * (deg + 1)] = double(-0.5);
                c[2 * (deg + 1) + 2] = double(1.5);
            }
            if(deg > 2)
            {
                double factor_1, factor_2;
                for(int n = 3; n != deg + 1; n++)
                {
                    factor_1 = double(2 * n - 1) / double(n);
                    factor_2 = -double(n - 1) / double(n);
                    for(int i = 0; i != n; i++)
                    {
                        c[n * (deg + 1) + 1 + i] = factor_1 * c[(n - 1) * (deg + 1) + i];
                        c[n * (deg + 1) + i] += factor_2 * c[(n - 2) * (deg + 1) + i];
                    }
                }
            }
        }

        /*Construct the Legendre polynomials using polynomial basis.
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
        static void fillLegendreCoefficients(float* c, int deg)
        {
            std::fill(c, &c[(deg + 1) * (deg + 1)], float(0.0));
            c[0] = float(1.0);
            if(deg > 0)
            {
                c[deg + 1 + 1] = float(1.0);
            }
            if(deg > 1)
            {
                c[2 * (deg + 1)] = float(-0.5);
                c[2 * (deg + 1) + 2] = float(1.5);
            }
            if(deg > 2)
            {
                float factor_1, factor_2;
                for(int n = 3; n != deg + 1; n++)
                {
                    factor_1 = float(2 * n - 1) / float(n);
                    factor_2 = -float(n - 1) / float(n);
                    for(int i = 0; i != n; i++)
                    {
                        c[n * (deg + 1) + 1 + i] = factor_1 * c[(n - 1) * (deg + 1) + i];
                        c[n * (deg + 1) + i] += factor_2 * c[(n - 2) * (deg + 1) + i];
                    }
                }
            }
        }

        /**Values of the Legendre polynomials at specific time point.
         *
         * This function needs to be protected by mutex due to filling array of powers.
         * Alternative implementation is to create local array of powers.
         */
        void valuesAt(double t, double* array) const override
        {
            if(t < start)
            {
                t = start;
            }
            if(t > end)
            {
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
                        += legendreCoefficientsD[i * (polynomialDegree + 1) + n] * xn[n];
                }
            }
            delete[] xn;
        }

        /**Values of the Legendre polynomials at specific time point.
         *
         */
        void valuesAt(double t, float* array) const override
        {
            if(t < start)
            {
                t = start;
            }
            if(t > end)
            {
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
                        += legendreCoefficientsF[i * (polynomialDegree + 1) + n] * xn[n];
                }
            }
            delete[] xn;
        }

        /**Into the array insert the polynomial basis of Legendre polynomial of certain
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
            std::copy(&legendreCoefficientsF[deg * (polynomialDegree + 1)],
                      &legendreCoefficientsF[deg * (polynomialDegree + 1) + deg + 1], array);
        }

        /**Into the array insert the polynomial basis of Legendre polynomial of certain
         * polynomialDegree.
         *
         * array must be prealocated to the size deg+1
         * deg must be between 0 and polynomialDegree including
         */
        void polynomialValues(uint32_t deg, double* array)
        {
            if( deg > polynomialDegree)
            {
                io::throwerr("Degree %d must be in [0, %d]!", deg, polynomialDegree);
            }
            std::copy(&legendreCoefficientsD[deg * (polynomialDegree + 1)],
                      &legendreCoefficientsD[deg * (polynomialDegree + 1) + deg + 1], array);
        }

        void printLegendrePolynomial(uint32_t deg)
        {
            if(deg > polynomialDegree)
            {
                io::throwerr("Degree %d must be in [0, %d]!", deg, polynomialDegree);
            }
            std::cout << io::xprintf("L_%d(x) = ", deg);
            if(legendreCoefficientsD[deg * (polynomialDegree + 1)] != 0)
            {
                std::cout << io::xprintf(" %.1f",
                                         legendreCoefficientsD[deg * (polynomialDegree + 1)]);
            }
            for(uint32_t i = 1; i < deg + 1; i++)
            {
                if(legendreCoefficientsD[deg * (polynomialDegree + 1) + i] != 0)
                {
                    double c = legendreCoefficientsD[deg * (polynomialDegree + 1) + i];
                    if(c < 0)
                    {
                        std::cout << io::xprintf(
                            " %.1fx^%d", legendreCoefficientsD[deg * (polynomialDegree + 1) + i],
                            i);
                    } else
                    {
                        std::cout << io::xprintf(
                            " +%.1fx^%d", legendreCoefficientsD[deg * (polynomialDegree + 1) + i],
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
        /**Function that transforms the value t on the interval [start, end] to the value t' on the
         * interval [-1,1] that is support of Legendre polynomials.
         *
         */
        double transformationSlope;
        double transformationIntercept;
        uint32_t polynomialDegree;
        uint32_t startReportingDegree;

        double* legendreCoefficientsD;
        float* legendreCoefficientsF;
        double* xnD;
        float* xnF;
        mutable std::mutex powerProtectionMutex;
    };

} // namespace util
} // namespace CTL
