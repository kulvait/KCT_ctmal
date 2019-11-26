#pragma once

#include "stringFormatter.h"

#ifdef DEBUG

#include "matplotlibcpp.h"
namespace plt = matplotlibcpp;

#endif

namespace CTL {
namespace util {
    /**Interface that will be implemented to contain actual vector function supported on the
     * interval [start, end]
     *
     * Might be used to provide values of the base function.
     */
    class VectorFunctionI
    {
    public:
        /**Initializes support interval.
         *
         *Dimension is the size of the vector represented by the VectorFunctionI derivative class.
         */
        VectorFunctionI(uint32_t dimension, double start, double end)
        {
            if(start < end)
            {
                this->start = start;
                this->end = end;
            } else
            {
                io::throwerr("VectorFunctionI:You have supplied the range [start, end] = [%f, %f] "
                             "but the value of "
                             "start needs to be less then the value of end.",
                             start, end);
            }
            if(dimension < 1)
            {
                io::throwerr("Dimension must be one or more but you have provided %d", dimension);
            }
            this->dimension = dimension;
        }

        /**Dimension of the vector.
         *
         */
        uint32_t getDimension() { return dimension; }

        /**Provides values of the function at particular point.
         *
         *Array vals is not allocated by function itself but it must be preallocated.
         *The number of elements that needs to be allocated can be obtained by calling
         *VectorFunctionI::getDimension() and subtracting startReportingDegree
         *
         * @param t Time that represents the point of the range [start, end] or possibly is outside
         *this range.
         * @param vals Pointer to the first element of the array to fill with dimension values.
         * @param startReportingDegree Index of the first function to report.
         */
        virtual void valuesAt(double t, double* vals) const = 0;

        /**Provides values of the function at particular point.
         *
         *Array vals is not allocated by function itself but it must be preallocated.
         *The number of elements that needs to be allocated can be obtained by calling
         *VectorFunctionI::getDimension() and subtracting startReportingDegree
         *
         * @param t Time that represents the point of the range [start, end] or possibly is outside
         *this range.
         * @param vals Pointer to the first element of the array to fill with dimension values.
         * @param startReportingDegree Index of the first function to report.
         */
        virtual void valuesAt(double t, float* vals) const = 0;

        /**Start of support.
         */
        double getStart() { return start; }

        /**Start of support.
         */
        double getEnd() { return end; }

#ifdef DEBUG
        virtual void plotFunctions(uint32_t granularity = 100,
                                   std::shared_ptr<std::vector<std::string>> names = nullptr)
        {
            if(granularity < 2)
            {
                io::throwerr("It is not possible to plot functions with granularity less than 2");
            }
            std::vector<double> taxis;
            std::vector<std::vector<double>> values;
            for(uint32_t j = 0; j != dimension; j++)
            {
                values.push_back(std::vector<double>());
            }
            float* val = new float[dimension];
            double time = start;
            double increment = (end - start) / double(granularity - 1);
            for(uint32_t i = 0; i != granularity; i++)
            {
                taxis.push_back(time);
                valuesAt(time, val);
                for(uint32_t j = 0; j != dimension; j++)
                {
                    values[j].push_back(val[j]);
                }
                time += increment;
            }
            for(uint32_t j = 0; j != dimension; j++)
            {
                if(names == nullptr)
                {
                    plt::named_plot(io::xprintf("Function %d", j), taxis, values[j]);
                } else
                {
                    plt::named_plot((*names)[j], taxis, values[j]);
                }
            }
            plt::legend();
            plt::show();
            delete[] val;
        }
#endif

    protected:
        double start; ///< Start of interval of the support.
        double end; ///< End of interval of the support.
        uint32_t dimension; ///< Dimension of the output.
    };

} // namespace util
} // namespace CTL
