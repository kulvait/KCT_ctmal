#pragma once

namespace CTL {
namespace util {
    /**Interface that will be implemented to contain actual scalar function supported on the
     * interval [start, end]
     *
     * Might be used to provide values of the base function.
     */
    class ScalarFunctionI
    {
    public:
        /**Initializes support interval.
         *
         */
        ScalarFunctionI(double start, double end)
        {
            if(start < end)
            {
                this->start = start;
                this->end = end;
            } else
            {
                io::throwerr("ScalarFunctionI: You have supplied the range [start, end] = [%f, %f] but the value of "
                             "start needs to be less then the value of end.",
                             start, end);
            }
        }

        /**Provides value of the function at particular point.
         *
         */
        virtual double valueAt(double t) const = 0;

    protected:
        double start; ///< Start of interval of the support.
        double end; ///< End of interval of the support.
    };

} // namespace util
} // namespace CTL
