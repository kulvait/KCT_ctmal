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

namespace KCT::util {
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
                 uint32_t startReportDegree = 0,
                 double* valuesBeforeStart = nullptr,
                 double* valuesAfterEnd = nullptr)
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
        if(valuesBeforeStart != nullptr)
        {
            this->valuesBeforeStart = new double[baseSize];
            std::copy(valuesBeforeStart, valuesBeforeStart + baseSize, this->valuesBeforeStart);
        } else
        {
            this->valuesBeforeStart = nullptr;
        }
        if(valuesAfterEnd != nullptr)
        {
            this->valuesAfterEnd = new double[baseSize];
            std::copy(valuesAfterEnd, valuesAfterEnd + baseSize, this->valuesAfterEnd);
        } else
        {
            this->valuesAfterEnd = nullptr;
        }

        // Now precompute the values of legendre polynomials
    } ///< Inits the function

    StepFunction(std::string inputFunctionsFile,
                 uint32_t basisSize,
                 double start,
                 double end,
                 double timeShift = 0.0,
                 uint32_t startReportDegree = 0,
                 double* valuesBeforeStart = nullptr,
                 double* valuesAfterEnd = nullptr)
        : VectorFunctionI(basisSize - startReportDegree, start, end)
        , baseSize(basisSize)
        , startReportDegree(startReportDegree)
        , timeShift(timeShift)
    {
        std::string ERR;
        if(startReportDegree < 0)
        {
            ERR = io::xprintf(
                "Variable startReportDegree must be non negative but supplied was %d.",
                startReportDegree);
            KCTERR(ERR);
        }
        if(baseSize - startReportDegree < 0)
        {
            ERR = io::xprintf("Base size must be greater then startReportDegree and not "
                              "%d as supplied.",
                              baseSize);
            KCTERR(ERR);
        }
        io::DenFileInfo bfi(inputFunctionsFile);
        if(basisSize > bfi.dimz())
        {
            ERR = io::xprintf(
                "Basis size %d is greater than the number of functions in file %s, %d!", basisSize,
                inputFunctionsFile.c_str(), bfi.dimz());
            KCTERR(ERR);
        }
        this->valuesPerFunction = bfi.dimx();
        // The following needs to be defined after previous line
        transformationSlope = (double(valuesPerFunction) - 1.0) / (end - start);
        transformationIntercept = -start * (double(valuesPerFunction) - 1.0) / (end - start);
        this->valuesD = new double[baseSize * valuesPerFunction];
        this->valuesF = new float[baseSize * valuesPerFunction];
        // Fill this array only by values without offset.
        if(bfi.getElementType() == io::DenSupportedType::FLOAT64)
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
        } else if(bfi.getElementType() == io::DenSupportedType::FLOAT32)
        {

            std::shared_ptr<io::Frame2DReaderI<float>> pr
                = std::make_shared<io::DenFrame2DReader<float>>(inputFunctionsFile);
            std::shared_ptr<io::Frame2DI<float>> f;
            for(uint32_t i = 0; i != baseSize; i++)
            {
                f = pr->readFrame(i);
                for(uint32_t j = 0; j != valuesPerFunction; j++)
                {
                    valuesD[i * valuesPerFunction + j] = f->get(j, 0);
                    valuesF[i * valuesPerFunction + j] = f->get(j, 0);
                }
            }
        } else
        {
            ERR = io::xprintf("Basis should be encoded in double or float file!");
            KCTERR(ERR);
        }
        if(valuesBeforeStart != nullptr)
        {
            this->valuesBeforeStart = new double[baseSize];
            std::copy(valuesBeforeStart, valuesBeforeStart + baseSize, this->valuesBeforeStart);
        } else
        {
            this->valuesBeforeStart = nullptr;
        }
        if(valuesAfterEnd != nullptr)
        {
            this->valuesAfterEnd = new double[baseSize];
            std::copy(valuesAfterEnd, valuesAfterEnd + baseSize, this->valuesAfterEnd);
        } else
        {
            this->valuesAfterEnd = nullptr;
        }
    } ///< Inits the function

    ~StepFunction()
    {
        delete[] valuesD;
        delete[] valuesF;
        if(valuesBeforeStart != nullptr)
        {
            delete[] valuesBeforeStart;
        }
        if(valuesAfterEnd != nullptr)
        {
            delete[] valuesAfterEnd;
        }
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

    friend void swap(StepFunction& x, StepFunction& y) // nothrow
    {
        using std::swap;
        swap(static_cast<VectorFunctionI&>(x), static_cast<VectorFunctionI&>(y));
        swap(x.baseSize, y.baseSize);
        swap(x.valuesPerFunction, y.valuesPerFunction);
        swap(x.timeShift, y.timeShift);
        swap(x.startReportDegree, y.startReportDegree);
        swap(x.valuesBeforeStart, y.valuesBeforeStart);
        swap(x.valuesAfterEnd, y.valuesAfterEnd);
    }

    StepFunction& operator=(StepFunction other)
    {
        swap(*this, other);
        return *this;
    } // Assignment

    StepFunction(StepFunction&& other) noexcept
        : StepFunction(nullptr,
                       0,
                       0,
                       0.0,
                       0.0) // initialize with 0 breakpoints just to be able to destruct it
    {
        swap(*this, other);
    }

    /**Values of the function at specific time point.
     *
     * This function needs to be protected by mutex due to filling array of powers.
     * Alternative implementation is to create local array of powers.
     */
    void valuesAt(double t, double* array) const override
    {
        if(t < start && valuesBeforeStart != nullptr)
        {
            for(uint32_t i = startReportDegree; i < baseSize; i++)
            {
                array[i - startReportDegree] = valuesBeforeStart[i];
            }
            return;
        }
        if(t > end && valuesAfterEnd != nullptr)
        {
            for(uint32_t i = startReportDegree; i < baseSize; i++)
            {
                array[i - startReportDegree] = valuesAfterEnd[i];
            }
            return;
        }
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
        if(t < start && valuesBeforeStart != nullptr)
        {
            for(uint32_t i = startReportDegree; i < baseSize; i++)
            {
                array[i - startReportDegree] = valuesBeforeStart[i];
            }
            return;
        }
        if(t > end && valuesAfterEnd != nullptr)
        {
            for(uint32_t i = startReportDegree; i < baseSize; i++)
            {
                array[i - startReportDegree] = valuesAfterEnd[i];
            }
            return;
        }
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
    /**Function that transforms the value t on the interval [start, end] to the value t' on
     * the interval [0,valuesPerFunction-1] to quickly evaluate function at any t.
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
    double* valuesBeforeStart;
    double* valuesAfterEnd;
};

} // namespace KCT::util
