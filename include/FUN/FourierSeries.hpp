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
     *written to the array. Polynomials start from the degree startReportDegree and ends by the
     *polynomial of the degree.
     *degree ... degree of polynomial to report last in the resulting vector
     * startReportDegree ... The degree of polynomial to report first in the resulting vector.
     */
    class FourierSeries : public VectorFunctionI
    {
    public:
        FourierSeries(int degree, double start, double end, int startReportDegree = 0)
            : VectorFunctionI(degree - startReportDegree, start, end)
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
                io::throwerr("Degree must be greater then startReportDegree and not "
                             "%d as supplied.",
                             degree);
            }
            this->degree = degree;
            this->startReportDegree = startReportDegree;
        } ///< Inits the function

        ~FourierSeries() {}

        FourierSeries(const FourierSeries& other)
            : FourierSeries(other.degree, other.start, other.end, other.startReportDegree)
        {

        } // Copy constructor

        FourierSeries& operator=(FourierSeries other)
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
                }
            }
            return *this;
        } // Assignment

#ifdef DEBUG
        void plotFunctions(uint32_t granularity = 100,
                           std::shared_ptr<std::vector<std::string>> names = nullptr) override
        {
            if(names == nullptr)
            {
                names = std::make_shared<std::vector<std::string>>();
                for(int i = startReportDegree; i <= degree; ++i)
                {
                    names->push_back(io::xprintf("Fourier %d", i));
                }
            }
            VectorFunctionI::plotFunctions(granularity, names);
        }
#endif

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
            double x = transformToSupport(t);
            int n;
            for(int i = startReportDegree; i < degree; i++)
            {
                n = (i+1) / 2;
                if(i == 0)
                {
                    array[i - startReportDegree] = 1.0;
                } else if(i % 2 == 1)
                {
                    array[i - startReportDegree] = std::cos(x * piF * n );
                } else
                {
                    array[i - startReportDegree] = std::sin(x * piF * n );
                }
            }
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
            double x = transformToSupport(t);
            int n;
            for(int i = startReportDegree; i < degree; i++)
            {
                n = (i+1) / 2;
                if(i == 0)
                {
                    array[i - startReportDegree] = 1.0;
                } else if(i % 2 == 1)
                {
                    array[i - startReportDegree] = std::cos(x * piD * n );
                } else
                {
                    array[i - startReportDegree] = std::sin(x * piD * n );
                }
            }
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
        const float piF = 3.141592653589793238462643383279502884;
        const float piD = 3.141592653589793238462643383279502884;

        int degree;
        int startReportDegree;
    };

} // namespace util
} // namespace CTL
