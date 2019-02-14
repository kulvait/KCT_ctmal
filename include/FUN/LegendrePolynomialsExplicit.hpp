#pragma once

#include <cmath>
#include <iostream>
#include <mutex>

#include "stringFormatter.h"
#include "VectorFunctionI.h"

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
    class LegendrePolynomialsExplicit : public VectorFunctionI
    {
    public:
        LegendrePolynomialsExplicit(int degree, double start, double end, int startReportDegree = 0)
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
            xnF = new float[degree + 1];
            xnD = new double[degree + 1];
            fillLegendreCoefficients(legendreCoefficientsD, degree);
            fillLegendreCoefficients(legendreCoefficientsF, degree);
        } ///< Inits the function

        ~LegendrePolynomialsExplicit()
        {
            delete[] legendreCoefficientsD;
            delete[] legendreCoefficientsF;
            delete[] xnF;
            delete[] xnD;
        }

        LegendrePolynomialsExplicit(const LegendrePolynomialsExplicit& other)
            : LegendrePolynomialsExplicit(
                  other.degree, other.start, other.end, other.startReportDegree)
        {

        } // Copy constructor

        LegendrePolynomialsExplicit& operator=(LegendrePolynomialsExplicit other)
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
                    xnF = new float[degree + 1];
                    xnD = new double[degree + 1];
                    legendreCoefficientsD = new double[(this->degree + 1) * (this->degree + 1)];
                    legendreCoefficientsF = new float[(this->degree + 1) * (this->degree + 1)];
                    fillLegendreCoefficients(legendreCoefficientsD, degree);
                    fillLegendreCoefficients(legendreCoefficientsF, degree);
                }
            }
            return *this;
        } // Assignment
        
#if DEBUG
        void plotFunctions(uint32_t granularity = 100, std::shared_ptr<std::vector<string>> names = nullptr) override
	{
		if(names == nullptr)
		{
			names = std::make_shared<std::vector<string>>();
			for(int i = startReportDegree; i <= degree; ++i)
			{
				names->push_back(io::xprintf("Legendre %d", i);
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
            // std::lock_guard<std::mutex> guard(powerProtectionMutex);//Big overhead
            double x = transformToSupport(t);
            double* xn = new double[degree + 1];
            xn[0] = 1.0;
            for(int n = 1; n < degree + 1; n++)
            {
                // xnD[n] = std::pow(x, double(n)); // Numerically more stable but slower
                xn[n] = xn[n - 1] * x;
            }
            std::fill(array, &array[degree - startReportDegree + 1], double(0.0));
            for(int i = startReportDegree; i < degree + 1; i++)
            {
                for(int n = 0; n != i + 1; n++)
                {
                    array[i - startReportDegree]
                        += legendreCoefficientsD[i * (degree + 1) + n] * xn[n];
                }
            }
            delete[] xn;
        }

        /**Values of the Legendre polynomials at specific time point.
         *
         */
        void valuesAt(double t, float* array) const override
        {
            // std::lock_guard<std::mutex> guard(powerProtectionMutex);//Big overhead
            double x = transformToSupport(t);
            double xk = 1.0;
            float* xn = new float[degree + 1];
            for(int n = 0; n < degree + 1; n++)
            {
                // xnF[n] = std::pow(x, float(n)); // Numerically more stable but slower
                xn[n] = float(
                    xk); // Power is double precission and particular value is single precision
                xk *= x;
            }
            std::fill(array, &array[degree - startReportDegree + 1], float(0.0));
            for(int i = startReportDegree; i < degree + 1; i++)
            {
                for(int n = 0; n != i + 1; n++)
                {
                    array[i - startReportDegree]
                        += legendreCoefficientsF[i * (degree + 1) + n] * xn[n];
                }
            }
            delete[] xn;
        }

        /**Into the array insert the polynomial basis of Legendre polynomial of certain degree.
         *
         * array must be prealocated to the size deg+1
         * deg must be between 0 and degree including
         */
        void polynomialValues(int deg, float* array)
        {
            if(deg < 0 || deg > degree)
            {
                io::throwerr("Degree %d must be in [0, %d]!", deg, degree);
            }
            std::copy(&legendreCoefficientsF[deg * (degree + 1)],
                      &legendreCoefficientsF[deg * (degree + 1) + deg + 1], array);
        }

        /**Into the array insert the polynomial basis of Legendre polynomial of certain degree.
         *
         * array must be prealocated to the size deg+1
         * deg must be between 0 and degree including
         */
        void polynomialValues(int deg, double* array)
        {
            if(deg < 0 || deg > degree)
            {
                io::throwerr("Degree %d must be in [0, %d]!", deg, degree);
            }
            std::copy(&legendreCoefficientsD[deg * (degree + 1)],
                      &legendreCoefficientsD[deg * (degree + 1) + deg + 1], array);
        }

        void printLegendrePolynomial(int deg)
        {
            if(deg < 0 || deg > degree)
            {
                io::throwerr("Degree %d must be in [0, %d]!", deg, degree);
            }
            std::cout << io::xprintf("L_%d(x) = ", deg);
            if(legendreCoefficientsD[deg * (degree + 1)] != 0)
            {
                std::cout << io::xprintf(" %.1f", legendreCoefficientsD[deg * (degree + 1)]);
            }
            for(int i = 1; i < deg + 1; i++)
            {
                if(legendreCoefficientsD[deg * (degree + 1) + i] != 0)
                {
                    double c = legendreCoefficientsD[deg * (degree + 1) + i];
                    if(c < 0)
                    {
                        std::cout << io::xprintf(" %.1fx^%d",
                                                 legendreCoefficientsD[deg * (degree + 1) + i], i);
                    } else
                    {
                        std::cout << io::xprintf(" +%.1fx^%d",
                                                 legendreCoefficientsD[deg * (degree + 1) + i], i);
                    }
                }
            }
            std::cout << std::endl;
        }

	//Transformation to the interval [start,end] => [-1,1]
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
        double* legendreCoefficientsD;
        float* legendreCoefficientsF;
        double* xnD;
        float* xnF;
        mutable std::mutex powerProtectionMutex;
    };

} // namespace util
} // namespace CTL
