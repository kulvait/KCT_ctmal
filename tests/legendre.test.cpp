#include "catch.hpp"
#include <cmath>
#include <plog/Log.h>

// Internal libs
#include "FUN/LegendrePolynomial.h"
#include "FUN/LegendrePolynomialsDerivatives.hpp"
#include "FUN/LegendrePolynomialsExplicit.hpp"
#include "stringFormatter.h"

using namespace CTL;

/// Testing correctness of legendre polynomials and its derivatives

TEST_CASE("CTL.util.LegendrePolynomialsExplicit.evaluate", "[legendre]")
{
    util::LegendrePolynomialsExplicit lp110(110, -1010.0, -20.3);
    for(int i = 0; i != 20; i++)
    {
        lp110.printLegendrePolynomial(i);
    }

    // Test support transformation
    double interval = -20.3 + 1010.0;
    double dt = double(interval / 19);
    double d = double(2.0 / 19);
    for(int i = 0; i != 20; i++)
    {

        double a = lp110.transformToSupport(-1010.0 + i * dt);
        //        LOGD << io::xprintf("Transformed value %f is %f.", -1010.0 + i * dt, a);
        REQUIRE(std::abs(a - double(d * i - 1.0)) < 0.000001);
    }

    double* values = new double[50];
    float* f = new float[10];
    // Check with https://en.wikipedia.org/wiki/Legendre_polynomials
    lp110.polynomialValues(10, values);
    REQUIRE(std::abs(values[0] - double(-63.0 / 256)) < 0.000001);
    REQUIRE(std::abs(values[1] - double(0.0)) < 0.000001);
    REQUIRE(std::abs(values[2] - double(3465.0 / 256)) < 0.000001);
    REQUIRE(std::abs(values[3] - double(0.0)) < 0.000001);
    REQUIRE(std::abs(values[4] - double(-30030.0 / 256)) < 0.000001);
    REQUIRE(std::abs(values[5] - double(0.0)) < 0.000001);
    REQUIRE(std::abs(values[6] - double(90090.0 / 256)) < 0.000001);
    REQUIRE(std::abs(values[7] - double(0.0)) < 0.000001);
    REQUIRE(std::abs(values[8] - double(-109395.0 / 256)) < 0.000001);
    REQUIRE(std::abs(values[9] - double(0.0)) < 0.000001);
    REQUIRE(std::abs(values[10] - double(46189.0 / 256)) < 0.000001);
    lp110.polynomialValues(9, f);
    REQUIRE(std::abs(f[0] - float(0.0 / 128)) < 0.000001);
    REQUIRE(std::abs(f[1] - float(315.0 / 128)) < 0.000001);
    REQUIRE(std::abs(f[2] - float(0.0 / 128)) < 0.000001);
    REQUIRE(std::abs(f[3] - float(-4620.0 / 128)) < 0.000001);
    REQUIRE(std::abs(f[4] - float(0.0 / 128)) < 0.000001);
    REQUIRE(std::abs(f[5] - float(18018.0 / 128)) < 0.000001);
    REQUIRE(std::abs(f[6] - float(0.0 / 128)) < 0.000001);
    REQUIRE(std::abs(f[7] - float(-25740.0 / 128)) < 0.000001);
    REQUIRE(std::abs(f[8] - float(0.0 / 128)) < 0.000001);
    REQUIRE(std::abs(f[9] - float(12155.0 / 128)) < 0.000001);
    lp110.polynomialValues(8, f);
    REQUIRE(std::abs(f[0] - float(35.0 / 128)) < 0.000001);
    REQUIRE(std::abs(f[1] - float(0.0 / 128)) < 0.000001);
    REQUIRE(std::abs(f[2] - float(-1260.0 / 128)) < 0.000001);
    REQUIRE(std::abs(f[3] - float(0.0 / 128)) < 0.000001);
    REQUIRE(std::abs(f[4] - float(6930.0 / 128)) < 0.000001);
    REQUIRE(std::abs(f[5] - float(0.0 / 128)) < 0.000001);
    REQUIRE(std::abs(f[6] - float(-12012.0 / 128)) < 0.00001);
    REQUIRE(std::abs(f[7] - float(0.0 / 128)) < 0.000001);
    REQUIRE(std::abs(f[8] - float(6435.0 / 128)) < 0.000001);
    lp110.polynomialValues(7, values);
    REQUIRE(std::abs(values[0] - double(0.0 / 256)) < 0.000001);
    REQUIRE(std::abs(values[1] - double(-35.0 / 16)) < 0.000001);
    REQUIRE(std::abs(values[2] - double(0.0 / 256)) < 0.000001);
    REQUIRE(std::abs(values[3] - double(315.0 / 16)) < 0.000001);
    REQUIRE(std::abs(values[4] - double(0.0 / 256)) < 0.000001);
    REQUIRE(std::abs(values[5] - double(-693.0 / 16)) < 0.000001);
    REQUIRE(std::abs(values[6] - double(0.0 / 256)) < 0.000001);
    REQUIRE(std::abs(values[7] - double(429.0 / 16)) < 0.000001);
    lp110.polynomialValues(6, f);
    REQUIRE(std::abs(f[0] - float(-5.0 / 16)) < 0.000001);
    REQUIRE(std::abs(f[1] - float(0.0 / 16)) < 0.000001);
    REQUIRE(std::abs(f[2] - float(105.0 / 16)) < 0.000001);
    REQUIRE(std::abs(f[3] - float(0.0 / 128)) < 0.000001);
    REQUIRE(std::abs(f[4] - float(-315.0 / 16)) < float(0.00001));
    REQUIRE(std::abs(f[5] - float(0.0 / 128)) < 0.000001);
    REQUIRE(std::abs(f[6] - float(231.0 / 16)) < 0.00001);
    lp110.polynomialValues(5, values);
    REQUIRE(std::abs(values[0] - double(0.0)) < 0.000001);
    REQUIRE(std::abs(values[1] - double(15.0 / 8)) < 0.000001);
    REQUIRE(std::abs(values[2] - double(0.0)) < 0.000001);
    REQUIRE(std::abs(values[3] - double(-70.0 / 8)) < 0.000001);
    REQUIRE(std::abs(values[4] - double(0.0)) < 0.000001);
    REQUIRE(std::abs(values[5] - double(63.0 / 8)) < 0.000001);
    lp110.polynomialValues(4, f);
    REQUIRE(std::abs(f[0] - float(3.0 / 8)) < 0.000001);
    REQUIRE(std::abs(f[1] - float(0.0)) < 0.000001);
    REQUIRE(std::abs(f[2] - float(-30.0 / 8)) < 0.000001);
    REQUIRE(std::abs(f[3] - float(0.0)) < 0.000001);
    REQUIRE(std::abs(f[4] - float(35.0 / 8)) < float(0.00001));
    lp110.polynomialValues(3, values);
    REQUIRE(std::abs(values[0] - double(0.0)) < 0.000001);
    REQUIRE(std::abs(values[1] - double(-3.0 / 2)) < 0.000001);
    REQUIRE(std::abs(values[2] - double(0.0)) < 0.000001);
    REQUIRE(std::abs(values[3] - double(5.0 / 2)) < 0.000001);
    lp110.polynomialValues(2, f);
    REQUIRE(std::abs(f[0] - float(-1.0 / 2)) < 0.000001);
    REQUIRE(std::abs(f[1] - float(0.0)) < 0.000001);
    REQUIRE(std::abs(f[2] - float(3.0 / 2)) < 0.000001);
    lp110.polynomialValues(1, values);
    REQUIRE(std::abs(values[0] - double(0.0)) < 0.000001);
    REQUIRE(std::abs(values[1] - double(1.0)) < 0.000001);
    lp110.polynomialValues(0, f);
    REQUIRE(std::abs(f[0] - float(1.0)) < 0.000001);
    // Shifted Legendre polynomial
    util::LegendrePolynomialsExplicit lp5(5, 0, 1.0);
    double t = 0.201;
    double a = lp5.transformToSupport(t);
    double targetSupportVal = t * 2 - 1; // transforming [0,1] -> [-1,1]
    REQUIRE((a - targetSupportVal) < 1e-30); // a should be equal to t*2-1
    // LOGD << io::xprintf("Evaluating at %f.", a);
    lp5.valuesAt(t, values, 0);
    double x[6];
    x[0] = 1.0;
    x[1] = 2 * t - 1.0;
    x[2] = 6 * t * t - 6 * t + 1.0;
    x[3] = 20 * t * t * t - 30 * t * t + 12 * t - 1;
    x[4] = 70 * t * t * t * t - 140 * t * t * t + 90 * t * t - 20 * t + 1;
    x[5] = 252 * t * t * t * t * t - 630 * t * t * t * t + 560 * t * t * t - 210 * t * t + 30 * t
        - 1;
    for(int i = 0; i != 6; i++)
    {
        //    LOGD << io::xprintf("L_%d(%f) = %f", i, t, values[i]);
        REQUIRE(std::abs(values[i] - x[i]) < 0.000001);
    }
    delete[] values;
    delete[] f;
}

/** Test that derivatives of first 5 derivatives of shifted Legendre polynomials agree.
 *
 */
TEST_CASE("CTL.util.LegendrePolynomialsDerivative.values", "[Legendre]")
{
    util::LegendrePolynomialsDerivatives ld5(5, 0.0, 1.0);
    double t = 0.201;
    double x[6];
    double* v1;
    v1 = new double[6];
    x[0] = 0.0;
    x[1] = 2.0;
    x[2] = 12 * t - 6.0;
    x[3] = 60 * t * t - 60 * t + 12.0;
    x[4] = 280 * t * t * t - 420 * t * t + 180 * t - 20.0;
    x[5] = 1260 * t * t * t * t - 2520 * t * t * t + 1680 * t * t - 420 * t + 30;
    ld5.valuesAt(t, v1, 0);
    for(int i = 0; i != 6; i++)
    {
        ld5.printLegendreDerivative(i);
        REQUIRE(std::abs(v1[0] - x[0]) < 1e-10);
        REQUIRE(std::abs(v1[1] - x[1]) < 1e-10);
        REQUIRE(std::abs(v1[2] - x[2]) < 1e-10);
        REQUIRE(std::abs(v1[3] - x[3]) < 1e-10);
        REQUIRE(std::abs(v1[4] - x[4]) < 1e-10);
        REQUIRE(std::abs(v1[5] - x[5]) < 1e-10);
    }
}

/**Test that derivative of polynomial L6 on support [-1,1] multiplyed by 2/(end-start) is the same
 *as the derivative of L6 on general support [start, end]
 */
TEST_CASE("CTL.util.LegendrePolynomialsDerivative.shifteval", "[Legendre]")
{
    double start = 0;
    double end = 1;
    double slope = 2.0 / (end - start);
    util::LegendrePolynomialsDerivatives ld(6, 0, 1.0);
    util::LegendrePolynomialsDerivatives lx(6, -1.0, 1.0); // Normal support
    double t = 0.201;
    double x = ld.transformToSupport(t);
    double v1, v2;
    ld.valuesAt(t, &v1, 6);
    lx.valuesAt(x, &v2, 6);
    REQUIRE(std::abs(v1 - v2 * slope) < 1e-10);
}

TEST_CASE("TEST:class LegendrePolynomial.", "[Polynomial]")
{
    util::LegendrePolynomial lp0(0, -1, 1);
    util::LegendrePolynomial lp1(1, -1, 1);
    util::LegendrePolynomial lp2(2, -1, 1);
    util::LegendrePolynomial lp3(3, -1, 1);
    util::LegendrePolynomial lp4(4, -1, 1);
    util::LegendrePolynomial lp5(5, -1, 1);

    double x[6];
    double t = 0.22;
    double val = lp5.valueAt(t);
    lp5.valuesAt(t, &x[0]);
    double difference = x[5] - val;
    REQUIRE(difference < 1e-30);
    // LOGD << io::xprintf("Value is %f and difference of the two methods is %f", value,
    // difference);
    REQUIRE(std::abs(x[5] - lp5.valueAt(t)) < 0.00001);
}
