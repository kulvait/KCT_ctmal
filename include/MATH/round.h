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
	float intsize = long(x)+1;
        return intsize+std::round(x-intsize);
    }
}

double roundLow(double x)
{
    if(x < double(0))
    {
        return std::round(x);
    } else
    {
	double intsize = long(x)+1;
        return intsize+std::round(x-intsize);
    }
}

long lroundLow(float x)
{
    if(x < float(0))
    {
        return std::lround(x);
    } else
    {
	float intsize = long(x)+1;
        return long(x)+1+std::lround(x-intsize);
    }
}

long lroundLow(double x)
{
    if(x < double(0))
    {
        return std::lround(x);
    } else
    {
	double intsize = long(x)+1;
        return long(x)+1+std::lround(x-intsize);
    }
}

} // namespace CTL::math
