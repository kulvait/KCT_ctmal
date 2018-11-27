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
    class StepFunction : public VectorFunctionI
    {
    public:
        StepFunction(double* values,
                     int baseSize,
                     int valuesPerFunction,
                     double start,
                     double end,
                     int startReportDegree = 0)
            : VectorFunctionI(baseSize - startReportDegree, start, end)
            , transformationSlope((double(valuesPerFunction)-1.0)/(end-start))
            , transformationIntercept(-start*(double(valuesPerFunction) - 1.0) / (end - start))
        {
            if(startReportDegree < 0)
            {
                io::throwerr("Variable startReportDegree must be non negative but supplied was %d.",
                             startReportDegree);
            }
            if(baseSize - startReportDegree < 0)
            {
                io::throwerr("Base size must be greater then startReportDegree and not "
                             "%d as supplied.",
                             baseSize);
            }
            this->baseSize = baseSize;
            this->valuesD = new double[baseSize * valuesPerFunction];
            this->valuesF = new float[baseSize * valuesPerFunction];
            this->valuesPerFunction = valuesPerFunction;
            for(int i = 0; i != baseSize; i++)
            {
                for(int j = 0; j != valuesPerFunction; j++)
                {
                    valuesD[i * valuesPerFunction + j] = values[i * valuesPerFunction + j];
                    valuesF[i * valuesPerFunction + j] = values[i * valuesPerFunction + j];
                }
            }
            this->startReportDegree = startReportDegree;

            // Now precompute the values of legendre polynomials
        } ///< Inits the function

        ~StepFunction()
        {
            delete[] valuesD;
            delete[] valuesF;
        }

        // Copy constructor
        StepFunction(const StepFunction& other)
            : StepFunction(other.valuesD,
                           other.baseSize,
                           other.valuesPerFunction,
                           other.start,
                           other.end,
                           other.startReportDegree)
        {
        }

        StepFunction& operator=(StepFunction other)
        {
            if(this != &other)
            {
                VectorFunctionI::operator=(other);
                this->startReportDegree = other.startReportDegree;
		this->start = other.start;
		this->end = other.end;
                if(this->baseSize != other.baseSize
                   || this->valuesPerFunction != other.valuesPerFunction)
                {
                    this->baseSize = other.baseSize;
                    this->valuesPerFunction = other.valuesPerFunction;
                    delete[] valuesF;
                    delete[] valuesD;
                    valuesF = new float[baseSize * valuesPerFunction];
                    valuesD = new double[baseSize * valuesPerFunction];
                }
                for(int i = 0; i != baseSize; i++)
                {
                    for(int j = 0; j != valuesPerFunction; j++)
                    {
                        valuesD[i * valuesPerFunction + j] = other.valuesD[i * valuesPerFunction + j];
                        valuesF[i * valuesPerFunction + j] = other.valuesF[i * valuesPerFunction + j];
                    }
                }
            }
            return *this;
        } // Assignment


        /**Values of the function at specific time point.
         *
         * This function needs to be protected by mutex due to filling array of powers.
         * Alternative implementation is to create local array of powers.
         */
        void valuesAt(double t, double* array) const override
        {
            // std::lock_guard<std::mutex> guard(powerProtectionMutex);//Big overhead
            int index = std::max(std::min(valuesPerFunction-1, int(std::lround(transformToIndex(t)))), 0);
            std::fill(array, &array[baseSize - startReportDegree], double(0.0));
            for(int i = startReportDegree; i < baseSize; i++)
            {
                    array[i - startReportDegree] = valuesD[valuesPerFunction*i+index];
            }
        }

        /**Values of the Legendre polynomials at specific time point.
         *
         */
        void valuesAt(double t, float* array) const override
        {
            // std::lock_guard<std::mutex> guard(powerProtectionMutex);//Big overhead
            int index = std::max(std::min(valuesPerFunction-1, int(std::lround(transformToIndex(t)))), 0);
            std::fill(array, &array[baseSize - startReportDegree], float(0.0));
            for(int i = startReportDegree; i < baseSize; i++)
            {
                    array[i - startReportDegree] = valuesD[valuesPerFunction*i+index];
            }
        }

        // Transformation to the interval [start,end] => [0,valuesPerFunction-1]
        double transformToIndex(double t) const
        {
            //return (t-start)*(double(valuesPerFunction)-1.0)/(end-start);
		return t*transformationSlope+transformationIntercept;
        }

    private:
        /**Function that transforms the value t on the interval [start, end] to the value t' on the
         * interval [-1,1] that is support of Legendre polynomials.
         *
         */
        double transformationSlope;
        double transformationIntercept;

	int valuesPerFunction;
        int baseSize;
        int startReportDegree;
        double start, end;
        double* valuesD;
        float* valuesF;
        mutable std::mutex powerProtectionMutex;
    };

} // namespace util
} // namespace CTL
