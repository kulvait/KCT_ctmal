#pragma once

#include <cmath>
#include <iostream>

#include "stringFormatter.h"
#include "FUN/LegendrePolynomialsExplicit.hpp"
#include "FUN/VectorFunctionI.h"

namespace CTL {
namespace util {
    /** Class for evaluation Legendre polynomials stretched from [-1,1] to a domain [start, end].
     *
     *Values at particular point of the domain are evaluated for Legendre polynomials. Values are
     *written to the array. Polynomials start from the degree startReportDegree and ends by the
     *polynomial of the degree.
     *degree ... degree of polynomial to report last in the resulting vector
     * startReportDegree ... The degree of polynomial to report first in the resulting vector.
     */
    class LegendrePolynomialsDerivatives : public VectorFunctionI
    {
    public:
        LegendrePolynomialsDerivatives(int degree,
                                       double start,
                                       double end,
                                       int startReportDegree = 0)
            : VectorFunctionI(degree + 1 - startReportDegree, start, end)
            , transformationSlope(double(2.0 / (end - start)))
            , transformationIntercept(-double((end + start) / (end - start)))
        {
            if(startReportDegree < 0)
            {
                io::throwerr("Variable startReportDegree must be non negative but supplied was %d.",
                             startReportDegree);
            }
            if(degree + 1 - startReportDegree < 0)
            {
                io::throwerr("Polynomial degree must be greater then startReportDegree and not "
                             "%d as supplied.",
                             degree);
            }
            this->degree = degree;
            this->startReportDegree = startReportDegree;
            // Now precompute the values of legendre polynomials
            legendreCoefficientsD = new double[(degree + 1) * (degree + 1)];
            legendreCoefficientsF = new float[(degree + 1) * (degree + 1)];
            derivativeCoefficientsD = new double[(degree + 1) * (degree + 1)];
            derivativeCoefficientsF = new float[(degree + 1) * (degree + 1)];
            xnF = new float[degree + 1];
            xnD = new double[degree + 1];
            LegendrePolynomialsExplicit::fillLegendreCoefficients(legendreCoefficientsD, degree);
            LegendrePolynomialsExplicit::fillLegendreCoefficients(legendreCoefficientsF, degree);
            fillLegendreDerivatives(legendreCoefficientsD, derivativeCoefficientsD, degree);
            fillLegendreDerivatives(legendreCoefficientsF, derivativeCoefficientsF, degree);
        } ///< Inits the function

        ~LegendrePolynomialsDerivatives()
        {
            delete[] legendreCoefficientsD;
            delete[] legendreCoefficientsF;
            delete[] derivativeCoefficientsD;
            delete[] derivativeCoefficientsF;
            delete[] xnF;
            delete[] xnD;
        }

        LegendrePolynomialsDerivatives(const LegendrePolynomialsDerivatives& other)
            : LegendrePolynomialsDerivatives(
                  other.degree, other.start, other.end, other.startReportDegree)
        {

        } // Copy constructor

        LegendrePolynomialsDerivatives& operator=(LegendrePolynomialsDerivatives other)
        {
            if(this != &other)
            {
                VectorFunctionI::operator=(other);
                this->startReportDegree = startReportDegree;
                if(this->degree != other.degree)
                {
                    this->degree = other.degree;
                    double interval = other.end - other.start;
                    this->transformationSlope = double(2.0 / interval);
                    this->transformationIntercept = -double((other.end + other.start) / interval);
                    delete[] legendreCoefficientsD;
                    delete[] legendreCoefficientsF;
                    delete[] xnF;
                    delete[] xnD;
                    delete[] derivativeCoefficientsD;
                    delete[] derivativeCoefficientsF;
                    legendreCoefficientsD = new double[(this->degree + 1) * (this->degree + 1)];
                    legendreCoefficientsF = new float[(this->degree + 1) * (this->degree + 1)];
                    derivativeCoefficientsD = new double[(degree + 1) * (degree + 1)];
                    derivativeCoefficientsF = new float[(degree + 1) * (degree + 1)];
                    xnF = new float[degree + 1];
                    xnD = new double[degree + 1];
                    LegendrePolynomialsExplicit::fillLegendreCoefficients(legendreCoefficientsD,
                                                                          degree);
                    LegendrePolynomialsExplicit::fillLegendreCoefficients(legendreCoefficientsF,
                                                                          degree);
                    fillLegendreDerivatives(legendreCoefficientsD, derivativeCoefficientsD, degree);
                    fillLegendreDerivatives(legendreCoefficientsF, derivativeCoefficientsF, degree);
                }
            }
            return *this;
        } // Assignment

        /**This function computes coeficients of polynomials of particular degrees to compute
         *derivative at particular point of interval [-1,1] that has the value relative to shifted
         *polynomial. To obtain values relative to [-1,1] is needed to divide these values by
         *transformationSlope.
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

        /**This function computes coeficients of polynomials of particular degrees to compute
         * derivative at particular point of interval [-1,1] that has the value relative to shifted
         * polynomial. To obtain values relative to [-1,1] is needed to divide these values by
         * transformationSlope.
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
            double x = transformToSupport(t);
            xnD[0] = 1.0;
            for(int n = 1; n != degree + 1; n++)
            {
                xnD[n] = xnD[n - 1] * x;
            }
            std::fill(array, &array[degree - startReportDegree + 1], double(0.0));
            for(int i = startReportDegree; i < degree + 1; i++)
            {
                for(int n = 0; n != i + 1; n++)
                {
                    array[i - startReportDegree]
                        += derivativeCoefficientsD[i * (degree + 1) + n] * xnD[n];
                }
            }
        }

        /**Values of the Legendre polynomials derivatives at specific time point.
         *
         */
        void valuesAt(double t, float* array) const override
        {
            float x = float(transformToSupport(t));
            xnF[0] = 1.0;
            for(int n = 1; n != degree + 1; n++)
            {
                xnF[n] = xnF[n - 1] * x;
            }
            std::fill(array, &array[degree - startReportDegree + 1], float(0.0));
            for(int i = startReportDegree; i < degree + 1; i++)
            {
                for(int n = 0; n != i + 1; n++)
                {
                    array[i - startReportDegree]
                        += derivativeCoefficientsF[i * (degree + 1) + n] * xnF[n];
                }
            }
        }

        /**Print formulas for derivatives of the polynomial of particular degree supported on
         * interval [-1,1].
         *
         */
        void printLegendreDerivative(int deg)
        {
            if(deg < 0 || deg > degree)
            {
                io::throwerr("Degree %d must be in [0, %d]!", deg, degree);
            }
            std::cout << io::xprintf("L_%d'(x) = ", deg);
            bool leadingSignPrinted = false;
            int e = 0;
            double c;
            while(e < deg + 1)
            {
                c = derivativeCoefficientsD[deg * (degree + 1) + e] / transformationSlope;
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

        void polynomialValues(int deg, float* array)
        {
            std::copy(&derivativeCoefficientsF[deg * (degree + 1)],
                      &derivativeCoefficientsF[deg * (degree + 1) + deg + 1], array);
        }

        void polynomialValues(int deg, double* array)
        {
            std::copy(&derivativeCoefficientsD[deg * (degree + 1)],
                      &derivativeCoefficientsD[deg * (degree + 1) + deg + 1], array);
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
        double transformationSlope;
        double transformationIntercept;

        int degree;
        int startReportDegree;
        double *legendreCoefficientsD, *derivativeCoefficientsD, *xnD;
        float *legendreCoefficientsF, *derivativeCoefficientsF, *xnF;
    };

} // namespace util
} // namespace CTL
