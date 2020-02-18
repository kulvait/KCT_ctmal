#pragma once

#include "mkl.h"
#include "rawop.h"
#include <cmath>

namespace CTL::math {
class SplineFitter
{
public:
    /** Constructs Spline Fitter.
     * @brief Creates the object for fitting splines. Follows the
     * https://software.intel.com/en-us/mkl-developer-reference-c-data-fitting-usage-examples
     * Deep copies all memory objects.
     *
     * @param numberOfPoints Number of intersection points to fit.
     * @param s_order See
     * https://software.intel.com/en-us/mkl-developer-reference-c-df-editppspline1d
     * @param s_type See
     * https://software.intel.com/en-us/mkl-developer-reference-c-df-editppspline1d
     */
    SplineFitter(MKL_INT numberOfBreakpoints,
                 MKL_INT s_order = DF_PP_CUBIC,
                 MKL_INT s_type = DF_PP_AKIMA)
        : numberOfBreakpoints(numberOfBreakpoints)
        , s_order(s_order)
        , s_type(s_type)
    {
        if(numberOfBreakpoints > 0)
        {
            if(s_order == DF_PP_LINEAR)
            {
                coeffs = new double[2 * (numberOfBreakpoints - 1)];
            } else if(s_order == DF_PP_QUADRATIC)
            {
                coeffs = new double[3 * (numberOfBreakpoints - 1)];
            } else if(s_order == DF_PP_CUBIC)
            {
                coeffs = new double[4 * (numberOfBreakpoints - 1)];
            }
            t_axe_internal = new double[numberOfBreakpoints];
            y_axe_internal = new double[numberOfBreakpoints];
        } else
        {
            coeffs = nullptr;
            t_axe_internal = nullptr;
            y_axe_internal = nullptr;
            bc_internal = nullptr;
        }
    }

    ~SplineFitter()
    {
        if(coeffs != nullptr)
        {
            delete[] coeffs;
            coeffs = nullptr;
        }
        if(t_axe_internal != nullptr)
        {
            delete[] t_axe_internal;
            t_axe_internal = nullptr;
        }
        if(y_axe_internal != nullptr)
        {
            delete[] y_axe_internal;
            y_axe_internal = nullptr;
        }
        if(bc_internal != nullptr)
        {
            delete[] bc_internal;
            bc_internal = nullptr;
        }
        if(task != nullptr)
        {
            dfDeleteTask(&task);
            task = nullptr;
        }
    }

    // Copy constructor
    SplineFitter(const SplineFitter& other)
        : SplineFitter(other.numberOfBreakpoints, other.s_order, other.s_type)
    {
        // If there is a task in other object already, lets copy it
        // using https://software.intel.com/en-us/forums/intel-math-kernel-library/topic/808247
        // Create new Data Fitting task
        // Set spline’s parameters and pointer to the precomputed spline coefficients as “scoeff”
        // parameter by df?EditPPSpline1D function (without spline reconstruction) Call
        // df?Interpolate1D function to run spline-based interpolation
        if(other.task != nullptr)
        {
            this->bc_type = other.bc_type;
            std::copy(other.t_axe_internal, other.t_axe_internal + numberOfBreakpoints,
                      t_axe_internal);
            std::copy(other.y_axe_internal, other.y_axe_internal + numberOfBreakpoints,
                      y_axe_internal);
            uint32_t coeffsSize;
            if(s_order == DF_PP_LINEAR)
            {
                coeffsSize = 2 * (numberOfBreakpoints - 1);
            } else if(s_order == DF_PP_QUADRATIC)
            {
                coeffsSize = 3 * (numberOfBreakpoints - 1);
            } else if(s_order == DF_PP_CUBIC)
            {
                coeffsSize = 4 * (numberOfBreakpoints - 1);
            }
            std::copy(other.coeffs, other.coeffs + coeffsSize, coeffs);
            if(other.bc_internal != nullptr)
            {
                bc_internal = new double[2];
                std::copy(other.bc_internal, other.bc_internal + 2, bc_internal);
            }
            MKL_INT status;
            // See https://software.intel.com/en-us/mkl-developer-reference-c-df-newtask1d
            status = dfdNewTask1D(&task, numberOfBreakpoints, t_axe_internal, DF_NO_HINT, 1,
                                  y_axe_internal, DF_NO_HINT);
            if(status != DF_STATUS_OK)
            {
                LOGE << io::xprintf("Can not create task, error %d! See mkl_df_defines.h", status);
            }
            // See https://software.intel.com/en-us/mkl-developer-reference-c-df-editppspline1d
            status = dfdEditPPSpline1D(task, s_order, s_type, bc_type, bc_internal, DF_NO_IC,
                                       nullptr, coeffs, DF_NO_HINT);
            if(status != DF_STATUS_OK)
            {
                LOGE << io::xprintf(
                    "Can not create spline interpolation, error %d! See mkl_df_defines.h", status);
            }
        }
    }

    friend void swap(SplineFitter& x, SplineFitter& y) // nothrow
    {
        using std::swap;
        swap(x.numberOfBreakpoints, y.numberOfBreakpoints);
        swap(x.s_order, y.s_order);
        swap(x.s_type, y.s_type);
        swap(x.coeffs, y.coeffs);
        swap(x.t_axe_internal, y.t_axe_internal);
        swap(x.y_axe_internal, y.y_axe_internal);
        swap(x.bc_internal, y.bc_internal);
        swap(x.task, y.task);
    }

    /**
     * Copy and swap see https://stackoverflow.com/questions/3279543/what-is-the-copy-and-swap-idiom
     *
     * @param other
     *
     * @return
     */
    SplineFitter& operator=(SplineFitter other)
    {
        swap(*this, other);
        return *this;
    }

    // move constructor
    SplineFitter(SplineFitter&& other) noexcept
        : SplineFitter(0) // initialize with 0 breakpoints just to be able to destruct it
    {
        swap(*this, other);
    }

    /**
     * This function do a spline fitting based on the data in the arrays t and y that are named
     * in https://software.intel.com/en-us/mkl-developer-reference-c-df-newtask1d t=x and y=y We
     * do a memory sanitization by copying memory into the object buffers since we don't know if
     * MKL internaly does deep copy to avoid possible memory leaks. Default parameters form invalid boundary condition.
     *
     * @param t
     * @param y
     * @param bc_type
     * @param bc In MKL documentation that is nullptr or one or two dimensional array based on
     * bc_type. For simplicity in this implementation we require to pass 2 dimensional double
     * array or nullptr even in the situation where MKL library require one dimensional array
     */
    void buildSpline(double* t,
                     double* y,
                     MKL_INT bc_type = DF_BC_1ST_LEFT_DER | DF_BC_1ST_RIGHT_DER,
                     double* bc = nullptr)
    {
        std::copy(t, t + numberOfBreakpoints, t_axe_internal);
        std::copy(y, y + numberOfBreakpoints, y_axe_internal);
        if(bc_internal != nullptr)
        {
            delete[] bc_internal;
        }
        if(bc != nullptr)
        {
            bc_internal = new double[2];
            std::copy(bc, bc + 2, bc_internal);
        } else
        {
            bc_internal = nullptr;
        }
        this->bc_type = bc_type;
        MKL_INT status;
        if(task != nullptr)
        {
            status = dfDeleteTask(&task);
            if(status != DF_STATUS_OK)
            {
                LOGE << "Can not dealocate previous task.";
            }
        }
        task = nullptr;
        // See https://software.intel.com/en-us/mkl-developer-reference-c-df-newtask1d
        status = dfdNewTask1D(&task, numberOfBreakpoints, t_axe_internal, DF_NO_HINT, 1,
                              y_axe_internal, DF_NO_HINT);
        if(status != DF_STATUS_OK)
        {
            LOGE << io::xprintf("Can not create task, error %d! See mkl_df_defines.h", status);
        }
        // See https://software.intel.com/en-us/mkl-developer-reference-c-df-editppspline1d
        status = dfdEditPPSpline1D(task, s_order, s_type, bc_type, bc, DF_NO_IC, nullptr, coeffs,
                                   DF_NO_HINT);
        if(status != DF_STATUS_OK)
        {
            LOGE << io::xprintf(
                "Can not create spline interpolation, error %d! See mkl_df_defines.h", status);
        }
        status = dfdConstruct1D(task, DF_PP_SPLINE, DF_METHOD_STD);
        if(status != DF_STATUS_OK)
        {
            if(status == DF_ERROR_BAD_NX)
            {
                LOGE << io::xprintf("Can not construct spline from task, DF_ERROR_BAD_NX ... "
                                    "Invalid number of breakpoints. For Akima splines, minimum 5 "
                                    "breakpoints is required.");
            } else if(status == DF_ERROR_BAD_IC)
            {
                LOGE << io::xprintf("Can not construct spline from task, DF_ERROR_BAD_IC ... Array "
                                    "of internal conditions for spline construction is not "
                                    "defined.");
            } else
            {
                LOGE << io::xprintf(
                    "Can not construct spline from task, error %d! See mkl_df_defines.h", status);
            }
        }
    }

    void interpolateAt(uint16_t npoints, double* t, double* vals, MKL_INT sitehint = DF_NO_HINT)
    {
        if(task != nullptr)
        {
            int dorder = 1;
            MKL_INT status
                = dfdInterpolate1D(task, DF_INTERP, DF_METHOD_PP, npoints, t, sitehint, 1, &dorder,
                                   DF_NO_APRIORI_INFO, vals, DF_MATRIX_STORAGE_ROWS, nullptr);
            if(status != DF_STATUS_OK)
            {
                LOGE << "Can not interpolate function.";
            }
        } else
        {
            LOGE << "Fit spline first";
        }
    }

private:
    MKL_INT numberOfBreakpoints;
    MKL_INT s_order;
    MKL_INT s_type;
    MKL_INT bc_type;
    double* coeffs = nullptr;
    double* t_axe_internal = nullptr;
    double* y_axe_internal = nullptr;
    double* bc_internal = nullptr;
    DFTaskPtr task = nullptr; // MKL object
};
} // namespace CTL::math
