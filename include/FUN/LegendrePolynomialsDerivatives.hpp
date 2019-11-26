#pragma once

#include <cmath>
#include <iostream>

#include "FUN/LegendrePolynomialsExplicit.hpp"
#include "FUN/VectorFunctionI.h"
#include "stringFormatter.h"

namespace CTL {
namespace util {
    /** Class for evaluation Legendre polynomials stretched from [-1,1] to a domain [start, end].
     *
     *Values at particular point of the domain are evaluated for Legendre polynomials. Values are
     *written to the array. Polynomials start from the polynomialDegree startReportDegree and ends
     *by the polynomial of the polynomialDegree. polynomialDegree ... polynomialDegree of polynomial
     *to report last in the resulting vector startReportDegree ... The polynomialDegree of
     *polynomial to report first in the resulting vector.
     */
    class LegendrePolynomialsDerivatives : public VectorFunctionI
    {
    public:
        LegendrePolynomialsDerivatives(uint32_t polynomialDegree,
                                       double start,
                                       double end,
                                       uint32_t startReportingDegree = 0)
            : VectorFunctionI(polynomialDegree + 1 - startReportingDegree, start, end)
            , polynomialDegree(polynomialDegree)
            , startReportingDegree(startReportingDegree)
            , transformationSlope(double(2.0 / (end - start)))
            , transformationIntercept(-double((end + start) / (end - start)))
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
            derivativeCoefficientsD = new double[(polynomialDegree + 1) * (polynomialDegree + 1)];
            derivativeCoefficientsF = new float[(polynomialDegree + 1) * (polynomialDegree + 1)];
            xnF = new float[polynomialDegree + 1];
            xnD = new double[polynomialDegree + 1];
            LegendrePolynomialsExplicit::fillLegendreCoefficients(legendreCoefficientsD,
                                                                  polynomialDegree);
            LegendrePolynomialsExplicit::fillLegendreCoefficients(legendreCoefficientsF,
                                                                  polynomialDegree);
            fillLegendreDerivatives(legendreCoefficientsD, derivativeCoefficientsD,
                                    polynomialDegree);
            fillLegendreDerivatives(legendreCoefficientsF, derivativeCoefficientsF,
                                    polynomialDegree);
        } ///< Inits the function

        ~LegendrePolynomialsDerivatives()
        {
            if(legendreCoefficientsD != nullptr)
            {
                delete[] legendreCoefficientsD;
            }
            if(legendreCoefficientsF != nullptr)
            {
                delete[] legendreCoefficientsF;
            }
            if(derivativeCoefficientsD != nullptr)
            {
                delete[] derivativeCoefficientsD;
            }
            if(derivativeCoefficientsF != nullptr)
            {
                delete[] derivativeCoefficientsF;
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
            derivativeCoefficientsD = nullptr;
            derivativeCoefficientsF = nullptr;
            xnF = nullptr;
            xnD = nullptr;
        }

        LegendrePolynomialsDerivatives(const LegendrePolynomialsDerivatives& other)
            : LegendrePolynomialsDerivatives(
                  other.polynomialDegree, other.start, other.end, other.startReportingDegree)
        {

        } // Copy constructor

        LegendrePolynomialsDerivatives& operator=(LegendrePolynomialsDerivatives other)
        {
            if(this != &other)
            {
                VectorFunctionI::operator=(other);
                this->startReportingDegree = other.startReportingDegree;
                double interval = other.end - other.start;
                this->transformationSlope = double(2.0 / interval);
                this->transformationIntercept = -double((other.end + other.start) / interval);
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
                    if(derivativeCoefficientsD != nullptr)
                    {
                        delete[] derivativeCoefficientsD;
                    }
                    if(derivativeCoefficientsF != nullptr)
                    {
                        delete[] derivativeCoefficientsF;
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
                    derivativeCoefficientsD = nullptr;
                    derivativeCoefficientsF = nullptr;
                    xnF = nullptr;
                    xnD = nullptr;
                    legendreCoefficientsD
                        = new double[(this->polynomialDegree + 1) * (this->polynomialDegree + 1)];
                    legendreCoefficientsF
                        = new float[(this->polynomialDegree + 1) * (this->polynomialDegree + 1)];
                    derivativeCoefficientsD
                        = new double[(polynomialDegree + 1) * (polynomialDegree + 1)];
                    derivativeCoefficientsF
                        = new float[(polynomialDegree + 1) * (polynomialDegree + 1)];
                    xnF = new float[polynomialDegree + 1];
                    xnD = new double[polynomialDegree + 1];
                    LegendrePolynomialsExplicit::fillLegendreCoefficients(legendreCoefficientsD,
                                                                          polynomialDegree);
                    LegendrePolynomialsExplicit::fillLegendreCoefficients(legendreCoefficientsF,
                                                                          polynomialDegree);
                    fillLegendreDerivatives(legendreCoefficientsD, derivativeCoefficientsD,
                                            polynomialDegree);
                    fillLegendreDerivatives(legendreCoefficientsF, derivativeCoefficientsF,
                                            polynomialDegree);
                }
            }
            return *this;
        } // Assignment

        /**This function computes coeficients of polynomials of particular polynomialDegrees to
         *compute derivative at particular point of interval [-1,1] that has the value relative to
         *shifted polynomial. To obtain values relative to [-1,1] is needed to divide these values
         *by transformationSlope.
         */
        void fillLegendreDerivatives(double* c, double* d, int deg)
        {
            std::fill(d, &d[(deg + 1) * (deg + 1)], double(0.0));
            for(int n = 0; n != deg + 1; n++)
                for(int i = 0; i != deg; i++)
                {
                    d[n * (deg + 1) + i]
                        = double(i + 1) * c[n * (deg + 1) + (i + 1)] * transformationSlope;
                }
        }

        /**This function computes coeficients of polynomials of particular polynomialDegrees to
         * compute derivative at particular point of interval [-1,1] that has the value relative to
         * shifted polynomial. To obtain values relative to [-1,1] is needed to divide these values
         * by transformationSlope.
         */
        void fillLegendreDerivatives(float* c, float* d, int deg)
        {
            std::fill(d, &d[(deg + 1) * (deg + 1)], float(0.0));
            for(int n = 0; n != deg + 1; n++)
                for(int i = 0; i != deg; i++)
                {
                    d[n * (deg + 1) + i]
                        = float(i + 1) * c[n * (deg + 1) + (i + 1)] * float(transformationSlope);
                }
        }

        /**Values of the Legendre polynomials derivatives at specific time point.
         *
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
            double x = transformToSupport(t);
            xnD[0] = 1.0;
            for(uint32_t n = 1; n != polynomialDegree + 1; n++)
            {
                xnD[n] = xnD[n - 1] * x;
            }
            std::fill(array, &array[polynomialDegree - startReportingDegree + 1], double(0.0));
            for(uint32_t i = startReportingDegree; i < polynomialDegree + 1; i++)
            {
                for(uint32_t n = 0; n != i + 1; n++)
                {
                    array[i - startReportingDegree]
                        += derivativeCoefficientsD[i * (polynomialDegree + 1) + n] * xnD[n];
                }
            }
        }

        /**Values of the Legendre polynomials derivatives at specific time point.
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
            float x = float(transformToSupport(t));
            xnF[0] = 1.0;
            for(uint32_t n = 1; n != polynomialDegree + 1; n++)
            {
                xnF[n] = xnF[n - 1] * x;
            }
            std::fill(array, &array[polynomialDegree - startReportingDegree + 1], float(0.0));
            for(uint32_t i = startReportingDegree; i < polynomialDegree + 1; i++)
            {
                for(uint32_t n = 0; n != i + 1; n++)
                {
                    array[i - startReportingDegree]
                        += derivativeCoefficientsF[i * (polynomialDegree + 1) + n] * xnF[n];
                }
            }
        }

        /**Print formulas for derivatives of the polynomial of particular polynomialDegree supported
         * on interval [-1,1].
         *
         */
        void printLegendreDerivative(uint32_t deg)
        {
            if(deg > polynomialDegree)
            {
                io::throwerr("Degree %d must be in [0, %d]!", deg, polynomialDegree);
            }
            std::cout << io::xprintf("L_%d'(x) = ", deg);
            bool leadingSignPrinted = false;
            uint32_t e = 0;
            double c;
            while(e < deg + 1)
            {
                c = derivativeCoefficientsD[deg * (polynomialDegree + 1) + e] / transformationSlope;
                if(!leadingSignPrinted)
                {
                    if(c != 0)
                    {
                        if(e == 0)
                        {
                            std::cout << io::xprintf(" %.1f", c);
                        } else if(e == 1)
                        {
                            std::cout << io::xprintf(" %.1fx", c);
                        } else
                        {
                            std::cout << io::xprintf(" %.1fx^%d", c, e);
                        }
                        leadingSignPrinted = true;
                    }
                } else
                {
                    if(c != 0)
                    {
                        if(c < 0)
                        {
                            std::cout << io::xprintf(" %.1fx^%d", c, e);
                        } else
                        {
                            std::cout << io::xprintf(" +%.1fx^%d", c, e);
                        }
                    }
                }
                e++;
            }
            if(!leadingSignPrinted)
            {
                std::cout << io::xprintf(" %.1f", 0.0);
            }
            std::cout << std::endl;
        }

        void polynomialValues(uint32_t deg, float* array)
        {
            std::copy(&derivativeCoefficientsF[deg * (polynomialDegree + 1)],
                      &derivativeCoefficientsF[deg * (polynomialDegree + 1) + deg + 1], array);
        }

        void polynomialValues(int deg, double* array)
        {
            std::copy(&derivativeCoefficientsD[deg * (polynomialDegree + 1)],
                      &derivativeCoefficientsD[deg * (polynomialDegree + 1) + deg + 1], array);
        }
        double transformToSupport(double t) const
        {
            return (transformationSlope * t) + transformationIntercept;
        }

    private:
        /**Function that transforms the value t on the interval [start, end] to the value t' on the
         * interval [-1,1] that is support of Legendre polynomials.
         *
         */
        uint32_t polynomialDegree;
        uint32_t startReportingDegree;
        double transformationSlope;
        double transformationIntercept;

        double *legendreCoefficientsD, *derivativeCoefficientsD, *xnD;
        float *legendreCoefficientsF, *derivativeCoefficientsF, *xnF;
    };

} // namespace util
} // namespace CTL
