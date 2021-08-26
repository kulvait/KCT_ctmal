#pragma once

#include <cmath>
#include <iostream>
#include <mutex>

#include "VectorFunctionI.h"
#include "stringFormatter.h"

namespace KCT {
namespace util {
    /** Class for evaluation of Fourier baisis stretched from [0,1] to a domain [start, end].
     *
     *Values at particular point of the domain are evaluated. Values are written to the array
     *consecutively starting from the function startReportDegree and ends by the number of Fourier
     *functions.
     *numberOfFunctions ... Number of functions that this class represents. First function is a
     *constant. Second function is sin(2*pi), third cos(2*pi), then sin(3*pi), cos(3*pi),
     *halfPeriodicFunctions... second functions are sin(pi), third cos(pi)
     */
    class FourierSeries : public VectorFunctionI
    {
    public:
        /**
         * @brief
         *
         * @param numberOfFunctions
         * @param start
         * @param end
         * @param halfPeriodicFunctions
         */
        FourierSeries(uint32_t numberOfFunctions,
                      double start,
                      double end,
                      uint32_t startReportingDegree = 0,
                      bool halfPeriodicFunctions = false,
                      bool constantOutsideInterval = true)
            : VectorFunctionI(numberOfFunctions - startReportingDegree, start, end)
            , transformationSlope(double(1.0 / (end - start)))
            , transformationIntercept(-start / (end - start))
            , numberOfFunctions(numberOfFunctions)
            , startReportingDegree(startReportingDegree)
            , halfPeriodicFunctions(halfPeriodicFunctions)
            , constantOutsideInterval(constantOutsideInterval)
        {
            if(startReportingDegree >= numberOfFunctions)
            {
                std::string msg
                    = io::xprintf("Parameter sartReportingDegree=%d but it should be less than "
                                  "numberOfFunctions=%d, correct!",
                                  startReportingDegree, numberOfFunctions);
                LOGD << msg;
                throw std::runtime_error(msg);
            }
        } ///< Inits the function

        ~FourierSeries() {}

        FourierSeries(const FourierSeries& other)
            : FourierSeries(other.numberOfFunctions,
                            other.start,
                            other.end,
                            other.startReportingDegree,
                            other.halfPeriodicFunctions)
        {

        } // Copy constructor

        FourierSeries& operator=(FourierSeries other)
        {
            if(this != &other)
            {
                VectorFunctionI::operator=(other);
                this->startReportingDegree = other.startReportingDegree;
                this->halfPeriodicFunctions = other.halfPeriodicFunctions;
                this->numberOfFunctions = other.numberOfFunctions;
                double interval = other.end - other.start;
                this->transformationSlope = double(1.0 / interval);
                this->transformationIntercept = double(-other.start / interval);
            }
            return *this;
        } // Assignment

#ifdef DEBUG
        void plotFunctions(uint32_t granularity = 100,
                           std::shared_ptr<std::vector<std::string>> names = nullptr,
                           std::string title = "") override
        {
            if(names == nullptr)
            {
                names = std::make_shared<std::vector<std::string>>();
                for(uint32_t i = startReportingDegree; i <= numberOfFunctions; ++i)
                {
                    names->push_back(io::xprintf("Funtion %d", i));
                }
            }
            if(title == "")
            {
                VectorFunctionI::plotFunctions(granularity, names, "Harmonic basis");
            } else
            {
                VectorFunctionI::plotFunctions(granularity, names, title);
            }
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
            double x = transformToSupport(t);
            int n;
            if(halfPeriodicFunctions)
            {
                for(uint32_t i = startReportingDegree; i < numberOfFunctions; i++)
                {
                    if(i == 0)
                    {
                        array[i - startReportingDegree] = 1.0;
                    } else if(i == 1)
                    {
                        array[i - startReportingDegree] = std::sin(x * piD);
                    } else if(i == 2)
                    {
                        array[i - startReportingDegree] = std::cos(x * piD);
                    } else
                    {
                        n = (i - 1) / 2;
                        if(i % 2 == 1)
                        {
                            array[i - startReportingDegree] = std::sin(x * 2 * piD * n);
                        } else
                        {
                            array[i - startReportingDegree] = std::cos(x * 2 * piD * n);
                        }
                    }
                }
            } else
            {
                for(uint32_t i = startReportingDegree; i < numberOfFunctions; i++)
                {
                    n = (i + 1) / 2;
                    if(i == 0)
                    {
                        array[i - startReportingDegree] = 1.0;
                    } else if(i % 2 == 1)
                    {
                        array[i - startReportingDegree] = std::sin(x * 2 * piD * n);
                    } else
                    {
                        array[i - startReportingDegree] = std::cos(x * 2 * piD * n);
                    }
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
            float x = transformToSupport(t);
            int n;
            if(halfPeriodicFunctions)
            {
                for(uint32_t i = startReportingDegree; i < numberOfFunctions; i++)
                {
                    if(i == 0)
                    {
                        array[i - startReportingDegree] = 1.0;
                    } else if(i == 1)
                    {
                        array[i - startReportingDegree] = std::sin(x * piD);
                    } else if(i == 2)
                    {
                        array[i - startReportingDegree] = std::cos(x * piD);
                    } else
                    {
                        n = (i - 1) / 2;
                        if(i % 2 == 1)
                        {
                            array[i - startReportingDegree] = std::sin(x * 2 * piD * n);
                        } else
                        {
                            array[i - startReportingDegree] = std::cos(x * 2 * piD * n);
                        }
                    }
                }
            } else
            {
                for(uint32_t i = startReportingDegree; i < numberOfFunctions; i++)
                {
                    n = (i + 1) / 2;
                    if(i == 0)
                    {
                        array[i - startReportingDegree] = 1.0;
                    } else if(i % 2 == 1)
                    {
                        array[i - startReportingDegree] = std::sin(x * 2 * piD * n);
                    } else
                    {
                        array[i - startReportingDegree] = std::cos(x * 2 * piD * n);
                    }
                }
            }
        }

        // Transformation to the interval [start,end] => [0,1]
        double transformToSupport(double t) const
        {
            return (transformationSlope * t) + transformationIntercept;
        }

    private:
        void valuesOutsideInterval(float* array) const
        {
            for(uint32_t i = startReportingDegree; i < numberOfFunctions; i++)
            {
                if(i % 2 == 0)
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
            for(uint32_t i = startReportingDegree; i < numberOfFunctions; i++)
            {
                if(i % 2 == 0)
                {
                    array[i - startReportingDegree] = 1.0;
                } else
                {
                    array[i - startReportingDegree] = 0.0;
                }
            }
        }
        /**Function that transforms the value [start, end] to [0,1].
         *
         */
        double transformationSlope;
        double transformationIntercept;
        const float piF = 3.141592653589793238462643383279502884;
        const float piD = 3.141592653589793238462643383279502884;

        uint32_t numberOfFunctions;
        uint32_t startReportingDegree;
        bool halfPeriodicFunctions;
        bool constantOutsideInterval;
    };

} // namespace util
} // namespace KCT
