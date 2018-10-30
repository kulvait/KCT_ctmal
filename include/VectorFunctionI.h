#pragma once

#include "stringFormatter.h"

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
         */
        VectorFunctionI(int dimension, double start, double end)
        {
            if(start < end)
            {
                this->start = start;
                this->end = end;
            } else
            {
                io::throwerr("You have supplied the range [start, end] = [%d, %d] but the value of "
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
        int getDimension() { return dimension; }

        /**Provides values of the function at particular point.
         *
         *Array vals is not allocated by function itself but it must be preallocated.
         *The number of elements that needs to be allocated can be obtained by calling
         *VectorFunctionI::getDimension()
         */
        virtual void valuesAt(double t, double* vals) const = 0;

        /**Provides values of the function at particular point.
         *
         *Array vals is not allocated by function itself but it must be preallocated.
         *The number of elements that needs to be allocated can be obtained by calling
         *VectorFunctionI::getDimension()
         */
        virtual void valuesAt(double t, float* vals) const = 0;

        /**Start of support.
         */
        double getStart() { return start; }

        /**Start of support.
         */
        double getEnd() { return end; }

    protected:
        double start; ///< Start of interval of the support.
        double end; ///< End of interval of the support.
        int dimension; ///< Dimension of the output.
    };

} // namespace util
} // namespace CTL