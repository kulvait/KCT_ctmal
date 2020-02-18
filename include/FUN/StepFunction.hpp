#pragma once

// External libraries
#include <cmath>
#include <iostream>
#include <mutex>

// Internal includes
#include "DEN/DenFrame2DReader.hpp"
#include "DEN/DenSupportedType.hpp"
#include "Frame2DReaderI.hpp"
#include "VectorFunctionI.h"
#include "stringFormatter.h"

namespace CTL {
namespace util {
    /** Class for evaluation of sampled functions on a domain [start, end].
     *
     *Values are written to the array. It is possible to omit first functions in the basis by
     *specifiing higher startReportDegree.
     *
     *@param[in] values Array with the function values
     *@param[in] baseSize number of the functions in the array
     *@param[in] valuesPerFunction number of sampling points
     *@param[in] start Start of interval.
     *@param[in] end End of interval.
     *@param[in] startReportDegree First function of the basis to report.
     */
    class StepFunction : public VectorFunctionI
    {
    public:
        StepFunction(double* values,
                     uint32_t baseSize,
                     uint32_t valuesPerFunction,
                     double start,
                     double end,
                     double timeShift = 0.0,
                     uint32_t startReportDegree = 0)
            : VectorFunctionI(baseSize - startReportDegree, start, end)
            , transformationSlope((double(valuesPerFunction) - 1.0) / (end - start))
            , transformationIntercept(-start * (double(valuesPerFunction) - 1.0) / (end - start))
            , timeShift(timeShift)
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
            for(uint32_t i = 0; i != baseSize; i++)
            {
                for(uint32_t j = 0; j != valuesPerFunction; j++)
                {
                    valuesD[i * valuesPerFunction + j] = values[i * valuesPerFunction + j];
                    valuesF[i * valuesPerFunction + j] = values[i * valuesPerFunction + j];
                }
            }
            this->startReportDegree = startReportDegree;

            // Now precompute the values of legendre polynomials
        } ///< Inits the function

        StepFunction(std::string inputFunctionsFile,
                     uint32_t basisSize,
                     double start,
                     double end,
                     double timeShift = 0.0,
                     uint32_t startReportDegree = 0)
            : VectorFunctionI(basisSize - startReportDegree, start, end)
            , baseSize(basisSize)
            , startReportDegree(startReportDegree)
            , timeShift(timeShift)
        {
            std::string err;
            if(startReportDegree < 0)
            {
                err = io::xprintf(
                    "Variable startReportDegree must be non negative but supplied was %d.",
                    startReportDegree);
                LOGE << err;
                throw std::runtime_error(err);
            }
            if(baseSize - startReportDegree < 0)
            {
                err = io::xprintf("Base size must be greater then startReportDegree and not "
                                  "%d as supplied.",
                                  baseSize);
                LOGE << err;
                throw std::runtime_error(err);
            }
            io::DenFileInfo bfi(inputFunctionsFile);
            if(basisSize > bfi.dimz())
            {
                err = io::xprintf(
                    "Basis size %d is greater than the number of functions in file %s, %d!",
                    basisSize, inputFunctionsFile.c_str(), bfi.dimz());
                LOGE << err;
                throw std::runtime_error(err);
            }
            this->valuesPerFunction = bfi.dimx();
            // The following needs to be defined after previous line
            transformationSlope = (double(valuesPerFunction) - 1.0) / (end - start);
            transformationIntercept = -start * (double(valuesPerFunction) - 1.0) / (end - start);
            this->valuesD = new double[baseSize * valuesPerFunction];
            this->valuesF = new float[baseSize * valuesPerFunction];
            // Fill this array only by values without offset.
            if(bfi.getDataType() == io::DenSupportedType::double_)
            {
                std::shared_ptr<io::Frame2DReaderI<double>> pr
                    = std::make_shared<io::DenFrame2DReader<double>>(inputFunctionsFile);
                std::shared_ptr<io::Frame2DI<double>> f;
                for(uint32_t i = 0; i != baseSize; i++)
                {
                    f = pr->readFrame(i);
                    for(uint32_t j = 0; j != valuesPerFunction; j++)
                    {
                        valuesD[i * valuesPerFunction + j] = f->get(j, 0);
                        valuesF[i * valuesPerFunction + j] = float(f->get(j, 0));
                    }
                }
            } else
            {
                err = io::xprintf("Basis should be encoded in double file!");
                LOGE << err;
                throw std::runtime_error(err);
            }
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
                           other.startReportDegree,
                           other.timeShift)
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
                for(uint32_t i = 0; i != baseSize; i++)
                {
                    for(uint32_t j = 0; j != valuesPerFunction; j++)
                    {
                        valuesD[i * valuesPerFunction + j]
                            = other.valuesD[i * valuesPerFunction + j];
                        valuesF[i * valuesPerFunction + j]
                            = other.valuesF[i * valuesPerFunction + j];
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
            t = t + timeShift;
            if(t < start)
            {
                t = start;
            }
            if(t > end)
            {
                t = end;
            }
            int index = std::max(
                std::min(int(valuesPerFunction) - 1, int(std::lround(transformToIndex(t)))), 0);
            for(uint32_t i = startReportDegree; i < baseSize; i++)
            {
                array[i - startReportDegree] = valuesD[valuesPerFunction * i + index];
            }
        }

        /**Values of the Legendre polynomials at specific time point.
         *
         */
        void valuesAt(double t, float* array) const override
        {
            t = t + timeShift;
            if(t < start)
            {
                t = start;
            }
            if(t > end)
            {
                t = end;
            }
            int index = std::max(
                std::min(int(valuesPerFunction) - 1, int(std::lround(transformToIndex(t)))), 0);
            for(uint32_t i = startReportDegree; i < baseSize; i++)
            {
                array[i - startReportDegree] = valuesF[valuesPerFunction * i + index];
            }
        }

        // Transformation to the interval [start,end] => [0,valuesPerFunction-1]
        double transformToIndex(double t) const
        {
            // return (t-start)*(double(valuesPerFunction)-1.0)/(end-start);
            return t * transformationSlope + transformationIntercept;
        }

    private:
        /**Function that transforms the value t on the interval [start, end] to the value t' on the
         * interval [0,valuesPerFunction-1] to quickly evaluate function at any t.
         *
         */
        double transformationSlope;
        double transformationIntercept;

        uint32_t baseSize;
        uint32_t startReportDegree;
        uint32_t valuesPerFunction;
        double* valuesD;
        float* valuesF;
        double timeShift;
    };

} // namespace util
} // namespace CTL
