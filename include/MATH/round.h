#pragma once

#include <cmath>

namespace CTL::math {

// Rounding in the way that
// roundLow(1.5) = 1
// roundLow(0.5) = 0
// roundLow(-0.5) = -1
float roundLow(float x)
{
    if(x < float(0))
    {
        return std::round(x);
    } else
    {
        return -std::round(-x);
    }
}

double roundLow(double x)
{
    if(x < double(0))
    {
        return std::round(x);
    } else
    {
        return -std::round(-x);
    }
}

long lroundLow(float x)
{
    if(x < float(0))
    {
        return std::lround(x);
    } else
    {
        return -std::lround(-x);
    }
}

long lroundLow(double x)
{
    if(x < double(0))
    {
        return std::lround(x);
    } else
    {
        return -std::lround(-x);
    }
}

} // namespace CTL::math
