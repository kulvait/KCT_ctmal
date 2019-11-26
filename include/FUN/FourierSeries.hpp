#pragma once

#include <cmath>
#include <iostream>
#include <mutex>

#include "VectorFunctionI.h"
#include "stringFormatter.h"

namespace CTL {
namespace util {
    /** Class for evaluation of Fourier baisis stretched from [0,1] to a domain [start, end].
     *
     *Values at particular point of the domain are evaluated. Values are written to the array
     *consecutively starting from the function startReportDegree and ends by the number of Fourier
     *functions.
     *degree ... Number of functions that this class represents. First function is a constant.
     *Second function is sin(2*pi), third cos(2*pi), then sin(3*pi), cos(3*pi),
     *halfPeriodicFunctions... second functions are sin(pi), third cos(pi)
     */
    class FourierSeries : public VectorFunctionI
    {
    public:
        /**
         * @brief
         *
         * @param degree
         * @param start
         * @param end
         * @param halfPeriodicFunctions
         */
        FourierSeries(int degree, double start, double end, bool halfPeriodicFunctions = false)
            : VectorFunctionI(degree, start, end)
            , transformationSlope(double(1.0 / (end - start)))
            , transformationIntercept(-start / (end - start))
        {
            this->degree = degree;
            this->startReportDegree = startReportDegree;
            this->halfPeriodicFunctions = halfPeriodicFunctions;
        } ///< Inits the function

        ~FourierSeries() {}

        FourierSeries(const FourierSeries& other)
            : FourierSeries(other.degree, other.start, other.end, other.halfPeriodicFunctions)
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
                           std::shared_ptr<std::vector<std::string>> names = nullptr,
                           uint32_t startReportingDegree = 0) override
        {
            if(names == nullptr)
            {
                names = std::make_shared<std::vector<std::string>>();
                for(int i = startReportDegree; i <= degree; ++i)
                {
                    names->push_back(io::xprintf("Fourier %d", i));
                }
            }
            VectorFunctionI::plotFunctions(granularity, names, startReportingDegree);
        }
#endif

        /**Values of the Legendre polynomials at specific time point.
         *
         * This function needs to be protected by mutex due to filling array of powers.
         * Alternative implementation is to create local array of powers.
         */
        void valuesAt(double t, double* array, uint32_t startReportingDegree) const override
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
            if(halfPeriodicFunctions)
            {
                for(int i = startReportDegree; i < degree; i++)
                {
                    if(i == 0)
                    {
                        array[i - startReportDegree] = 1.0;
                    } else if(i == 1)
                    {
                        array[i - startReportDegree] = std::sin(x * piD);
                    } else if(i == 2)
                    {
                        array[i - startReportDegree] = std::cos(x * piD);
                    } else
                    {
                        n = (i - 1) / 2;
                        if(i % 2 == 1)
                        {
                            array[i - startReportDegree] = std::sin(x * 2 * piD * n);
                        } else
                        {
                            array[i - startReportDegree] = std::cos(x * 2 * piD * n);
                        }
                    }
                }
            } else
            {
                for(int i = startReportDegree; i < degree; i++)
                {
                    n = (i + 1) / 2;
                    if(i == 0)
                    {
                        array[i - startReportDegree] = 1.0;
                    } else if(i % 2 == 1)
                    {
                        array[i - startReportDegree] = std::sin(x * 2 * piD * n);
                    } else
                    {
                        array[i - startReportDegree] = std::cos(x * 2 * piD * n);
                    }
                }
            }
        }

        /**Values of the Legendre polynomials at specific time point.
         *
         */
        void valuesAt(double t, float* array, uint32_t startReportingDegree) const override
        {
            if(t < start)
            {
                t = start;
            }
            if(t > end)
            {
                t = end;
            }
            float x = transformToSupport(t);
            int n;
            if(halfPeriodicFunctions)
            {
                for(int i = startReportDegree; i < degree; i++)
                {
                    if(i == 0)
                    {
                        array[i - startReportDegree] = 1.0;
                    } else if(i == 1)
                    {
                        array[i - startReportDegree] = std::sin(x * piD);
                    } else if(i == 2)
                    {
                        array[i - startReportDegree] = std::cos(x * piD);
                    } else
                    {
                        n = (i - 1) / 2;
                        if(i % 2 == 1)
                        {
                            array[i - startReportDegree] = std::sin(x * 2 * piD * n);
                        } else
                        {
                            array[i - startReportDegree] = std::cos(x * 2 * piD * n);
                        }
                    }
                }
            } else
            {
                for(int i = startReportDegree; i < degree; i++)
                {
                    n = (i + 1) / 2;
                    if(i == 0)
                    {
                        array[i - startReportDegree] = 1.0;
                    } else if(i % 2 == 1)
                    {
                        array[i - startReportDegree] = std::sin(x * 2 * piD * n);
                    } else
                    {
                        array[i - startReportDegree] = std::cos(x * 2 * piD * n);
                    }
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
        bool halfPeriodicFunctions;
    };

} // namespace util
} // namespace CTL
