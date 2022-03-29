#include "helpers.test.hpp"

std::array<double, 3> vecnorm(std::array<double, 3> v)
{
    std::array<double, 3> n;
    double na = std::sqrt(v[1] * v[1] + v[2] * v[2] + v[0] * v[0]);
    n[0] = v[0] / na;
    n[1] = v[1] / na;
    n[2] = v[2] / na;
    return n;
}

bool almostEqual(double x, double y, double tol)
{
    double dif = x - y;
    if(std::abs(dif) < tol)
    {
        return true;
    } else
    {
        return false;
    }
}
