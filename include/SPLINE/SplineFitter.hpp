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
    }

    ~SplineFitter()
    {
        if(coeffs != nullptr)
        {
            delete[] coeffs;
        }
        if(task != nullptr)
        {
            dfDeleteTask(&task);
        }
    }

    void buildSpline(double* t,
                     double* y,
                     MKL_INT bc_type = DF_BC_1ST_LEFT_DER | DF_BC_1ST_RIGHT_DER,
                     double* bc = nullptr)
    {
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
        status = dfdNewTask1D(&task, numberOfBreakpoints, t, DF_NO_HINT, 1, y, DF_NO_HINT);
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
    const MKL_INT numberOfBreakpoints;
    const MKL_INT s_order;
    const MKL_INT s_type;
    double* coeffs = nullptr;
    DFTaskPtr task = nullptr; // MKL object
};
} // namespace CTL::math
