#include "catch.hpp"
#include <cmath>
#include <plog/Log.h>

#define DEBUG

// Internal libs
#include "FUN/FourierSeries.hpp"
#include "matplotlibcpp.h"
#include "stringFormatter.h"
namespace plt = matplotlibcpp;

using namespace KCT;
TEST_CASE("KCT.util.FourierSeries", "[fourier]")
{
    util::FourierSeries ff(5, 1.0, 2.0, 0);
    ff.plotFunctions();
}
