#pragma once

// External libraries
#include <cmath>
#include <iostream>
#include <mutex>

// Internal includes
#include "DEN/DenFrame2DReader.hpp"
#include "DEN/DenSupportedType.hpp"
#include "Frame2DReaderI.hpp"
#include "SPLINE/SplineFitter.hpp"
#include "VectorFunctionI.h"
#include "stringFormatter.h"

namespace KCT {
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
     *@param[in] valuesBeforeStart Values to report before start of interval, if nullptr value at
     *start is used.
     *@param[in] startReportDegree Values to report after end of interval, if nullptr value ast end
     *is used.
     */
    class SplineInterpolatedFunction : public VectorFunctionI
    {
    public:
        SplineInterpolatedFunction(double* values,
                                   uint32_t baseSize,
                                   uint32_t valuesPerFunction,
                                   double start,
                                   double end,
                                   double timeShift = 0.0,
                                   uint32_t startReportDegree = 0,
                                   double* valuesBeforeStart = nullptr,
                                   double* valuesAfterEnd = nullptr)
            : VectorFunctionI(baseSize - startReportDegree, start, end)
            , baseSize(baseSize)
            , valuesPerFunction(valuesPerFunction)
            , timeShift(timeShift)
            , startReportDegree(startReportDegree)
        {
            if(baseSize == 0)
            {
                return;
            }
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
            // Time shift does not apply here but when evaluating functions
            defaultBoundaryConditions = new double[2];
            defaultBoundaryConditions[0] = 0.0;
            defaultBoundaryConditions[1] = 0.0;
            bc_type = DF_BC_1ST_LEFT_DER | DF_BC_1ST_RIGHT_DER;
            double* times = new double[valuesPerFunction];
            for(uint32_t j = 0; j != valuesPerFunction; j++)
            {
                times[j] = start + j * (end - start) / (double(valuesPerFunction) - 1.0);
            }
            this->valuesPerFunction = valuesPerFunction;
            std::shared_ptr<math::SplineFitter> sf;
            for(uint32_t i = 0; i != baseSize; i++)
            {
                sf = std::make_shared<math::SplineFitter>(valuesPerFunction);
                splineFitters.push_back(sf);
                splineFitters[i]->buildSpline(times, &values[i * valuesPerFunction], bc_type,
                                              defaultBoundaryConditions);
            }
            delete[] times;
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

        SplineInterpolatedFunction(std::string inputFunctionsFile,
                                   uint32_t basisSize,
                                   double start,
                                   double end,
                                   double timeShift = 0.0,
                                   uint32_t startReportDegree = 0,
                                   double* valuesBeforeStart = nullptr,
                                   double* valuesAfterEnd = nullptr)
            : VectorFunctionI(basisSize - startReportDegree, start, end)
            , baseSize(basisSize)
            , timeShift(timeShift)
            , startReportDegree(startReportDegree)
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
                KCTERR(err);
            }
            this->valuesPerFunction = bfi.dimx();
            // Fill this array only by values without offset.
            if(bfi.getElementType() == io::DenSupportedType::FLOAT64)
            {
                // Time shift does not apply here but when evaluating functions
                double* times = new double[valuesPerFunction];
                double* values = new double[valuesPerFunction];
                defaultBoundaryConditions = new double[2];
                defaultBoundaryConditions[0] = 0.0;
                defaultBoundaryConditions[1] = 0.0;
                bc_type = DF_BC_1ST_LEFT_DER | DF_BC_1ST_RIGHT_DER;
                for(uint32_t j = 0; j != valuesPerFunction; j++)
                {
                    times[j] = start + j * (end - start) / (double(valuesPerFunction) - 1.0);
                }
                std::shared_ptr<io::Frame2DReaderI<double>> pr
                    = std::make_shared<io::DenFrame2DReader<double>>(inputFunctionsFile);
                std::shared_ptr<io::Frame2DI<double>> f;
                std::shared_ptr<math::SplineFitter> sf;
                for(uint32_t i = 0; i != baseSize; i++)
                {
                    sf = std::make_shared<math::SplineFitter>(valuesPerFunction);
                    splineFitters.push_back(sf);
                    f = pr->readFrame(i);
                    for(uint32_t j = 0; j != valuesPerFunction; j++)
                    {
                        values[j] = f->get(j, 0);
                    }
                    splineFitters[i]->buildSpline(times, values, bc_type,
                                                  defaultBoundaryConditions);
                }
                delete[] times;
                delete[] values;
                if(valuesBeforeStart != nullptr)
                {
                    this->valuesBeforeStart = new double[baseSize];
                    std::copy(valuesBeforeStart, valuesBeforeStart + baseSize,
                              this->valuesBeforeStart);
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
            } else
            {
                err = io::xprintf("Basis should be encoded in double file!");
                LOGE << err;
                throw std::runtime_error(err);
            }
        } ///< Inits the function

        ~SplineInterpolatedFunction()
        {
            splineFitters.clear();
            if(defaultBoundaryConditions != nullptr)
            {
                delete[] defaultBoundaryConditions;
            }
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
        SplineInterpolatedFunction(const SplineInterpolatedFunction& other)
            : VectorFunctionI(other.baseSize - other.startReportDegree, other.start, other.end)
            , baseSize(other.baseSize)
            , valuesPerFunction(other.valuesPerFunction)
            , timeShift(other.timeShift)
            , startReportDegree(other.startReportDegree)
            , defaultBoundaryConditions(other.defaultBoundaryConditions)
            , bc_type(other.bc_type)
        {

            if(other.defaultBoundaryConditions != nullptr)
            {
                std::copy(other.defaultBoundaryConditions, other.defaultBoundaryConditions + 2,
                          defaultBoundaryConditions);
            }
            if(other.valuesBeforeStart != nullptr)
            {
                this->valuesBeforeStart = new double[baseSize];
                std::copy(other.valuesBeforeStart, other.valuesBeforeStart + baseSize,
                          valuesBeforeStart);
            }
            if(other.valuesAfterEnd != nullptr)
            {
                this->valuesAfterEnd = new double[baseSize];
                std::copy(other.valuesAfterEnd, other.valuesAfterEnd + baseSize, valuesAfterEnd);
            }
            std::shared_ptr<math::SplineFitter> sfa;
            for(std::shared_ptr<math::SplineFitter> const& sfb : other.splineFitters)
            {
                sfa = std::make_shared<math::SplineFitter>(*sfb); // Copy constructor
                this->splineFitters.push_back(sfa);
            }
        }

        friend void swap(SplineInterpolatedFunction& x, SplineInterpolatedFunction& y) // nothrow
        {
            using std::swap;
            swap(static_cast<VectorFunctionI&>(x), static_cast<VectorFunctionI&>(y));
            swap(x.baseSize, y.baseSize);
            swap(x.valuesPerFunction, y.valuesPerFunction);
            swap(x.timeShift, y.timeShift);
            swap(x.startReportDegree, y.startReportDegree);
            swap(x.splineFitters, y.splineFitters);
            swap(x.defaultBoundaryConditions, y.defaultBoundaryConditions);
            swap(x.bc_type, y.bc_type);
            swap(x.valuesBeforeStart, y.valuesBeforeStart);
            swap(x.valuesAfterEnd, y.valuesAfterEnd);
        }

        /**
         * Copy and swap see
         * https://stackoverflow.com/questions/3279543/what-is-the-copy-and-swap-idiom
         *
         * @param other
         *
         * @return
         */
        SplineInterpolatedFunction& operator=(SplineInterpolatedFunction other)
        {
            swap(*this, other);
            return *this;
        } // Assignment

        // move constructor
        SplineInterpolatedFunction(SplineInterpolatedFunction&& other) noexcept
            : SplineInterpolatedFunction(
                nullptr,
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
            double time = t;
            double v;
            for(uint32_t i = startReportDegree; i < baseSize; i++)
            {
                splineFitters[i]->interpolateAt(1, &time, &v);
                array[i - startReportDegree] = v;
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
            double v;
            for(uint32_t i = startReportDegree; i < baseSize; i++)
            {
                splineFitters[i]->interpolateAt(1, &t, &v);
                array[i - startReportDegree] = v;
            }
        }

    private:
        uint32_t baseSize;
        uint32_t valuesPerFunction;

        double timeShift;
        uint32_t startReportDegree;
        std::vector<std::shared_ptr<math::SplineFitter>> splineFitters;
        double* defaultBoundaryConditions = nullptr;
        MKL_INT bc_type = 0;
        double* valuesBeforeStart = nullptr;
        double* valuesAfterEnd = nullptr;
    };

} // namespace util
} // namespace KCT
