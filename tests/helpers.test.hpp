#pragma once

#include <iostream>

#include "MATRIX/utils.hpp"
#include "stringFormatter.h"

using namespace KCT;
using namespace KCT::matrix;

// Helper functions for ctmal testing
std::array<double, 3> vecnorm(std::array<double, 3> v);

bool almostEqual(double x, double y, double tol = 1e-10);

template <uint32_t N>
double normdiff(std::array<double, N> v, std::array<double, N> w)
{
    double nrmsq = 0.0;
    for(uint32_t i = 0; i != N; i++)
    {
        nrmsq += (v[i] - w[i]) * (v[i] - w[i]);
    }
    return std::sqrt(nrmsq);
}


template <int N>
bool almostEqual(std::array<double, N> x, std::array<double, N> y, double tol = 1e-10)
{
    std::array<double, N> minusx = multiplyVectorByConstant<N>(x, -1.0);
    std::array<double, N> vdf = vectorSum<N>(minusx, y);
    double n = vectorNorm<N>(vdf);
    if(n < tol)
    {
        return true;
    } else
    {
        return false;
    }
}

template <int N>
bool almostEqualRelative(std::array<double, N> x, std::array<double, N> y, double tol = 1e-10)
{
    std::array<double, N> minusx = multiplyVectorByConstant<N>(x, -1.0);
    std::array<double, N> vdf = vectorSum<N>(minusx, y);
    double basenorm = vectorNorm<N>(x);
    double n = vectorNorm<N>(vdf);
    if(n / basenorm < tol)
    {
        return true;
    } else
    {
        return false;
    }
}

template <int N>
std::string printVector(std::string name, std::array<double, N> x)
{
    std::string s = io::xprintf("%s: (", name.c_str());
    for(uint32_t i = 0; i < N - 1; i++)
    {
        s += io::xprintf("%f, ", x[i]);
    }
    s += io::xprintf("%f)", x[N - 1]);
    return s;
}
