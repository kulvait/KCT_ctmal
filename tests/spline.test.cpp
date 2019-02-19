#include "catch.hpp"
#include <cmath>
#include <plog/Log.h>

// Internal libs
#include "SPLINE/SplineFitter.hpp"
#include "stringFormatter.h"
#include "matplotlibcpp.h"
namespace plt = matplotlibcpp;

using namespace CTL;

/// Testing correctness of legendre polynomials and its derivatives

TEST_CASE("CTL.math.SplineFitter", "[spline]")
{
    int points = 6;
    math::SplineFitter fitter(points, DF_PP_CUBIC, DF_PP_AKIMA);
    double t[] = { 0.0, 1.0, 2.0, 3.0, 4.0, 5.0 };
    double y[] = { 0.0, 1.0, 0.0, 2.0, 0.0, 1.0 };
    MKL_INT bc_type = DF_BC_1ST_LEFT_DER | DF_BC_1ST_RIGHT_DER;
    double bc[] = { 0.0, 0.0 };
    fitter.buildSpline(t, y, bc_type, bc);

    // Now we produce the discretization
    uint32_t granularity = 100;
    std::vector<double> taxis;
    double* val = new double[granularity];
    double increment = 5.0 / (granularity - 1);
    for(uint32_t i = 0; i != granularity; i++)
    {
        taxis.push_back(i * increment);
    }
    fitter.interpolateAt(granularity, taxis.data(), val);
    std::vector<double> plotme;
    for(uint32_t i = 0; i != granularity; i++)
    {
        plotme.push_back(val[i]);
    }
    plt::plot(taxis, plotme);
    plt::show();
    delete[] val;
}

