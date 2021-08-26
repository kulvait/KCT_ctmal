#include "catch.hpp"
#include <cmath>
#include <plog/Log.h>

// Internal libs
#include "FUN/ChebyshevPolynomialsExplicit.hpp"
#include "stringFormatter.h"

using namespace KCT;

/// Testing correctness of legendre polynomials and its derivatives

TEST_CASE("KCT.util.ChebyshevPolynomialsExplicit.evaluate",
          "[polynomial][chebyshev][NOPRINT][NOVIZ]")
{
    // First define Chebyshev Polynomial of degree 110
    util::ChebyshevPolynomialsExplicit cp110(110, -1, 1, true, 0);
    //    for(int i = 0; i != 20; i++)
    //    {
    //        cp110.printChebyshevPolynomial(i);
    //    }

    double* values = new double[111];
    // Check with https://www.wolframalpha.com/input/?i=+chebyshevT%28110%2Cx%29
    cp110.polynomialValues(110, values);
    REQUIRE(std::abs(values[0] - -1.0) < 0.000001);
    REQUIRE(std::abs(values[1] - 0.0) < 0.000001);
    REQUIRE(std::abs(values[2] - 6050) < 0.000001);
    REQUIRE(std::abs(values[3] - 0.0) < 0.000001);
    REQUIRE(std::abs(values[4] - -6098400.0) < 0.000001);
    REQUIRE(std::abs(values[5] - 0.0) < 0.000001);
    REQUIRE(std::abs(values[6] - 2456435520.0) < 0.000001);
    REQUIRE(std::abs(values[7] - 0.0) < 0.000001);
    REQUIRE(std::abs(values[8] - -529186394880.0) < 0.000001);
    REQUIRE(std::abs(values[9] - 0.0) < 0.000001);
    REQUIRE(std::abs(values[10] - 70769860541952) < 0.000001);
    delete[] values;
}

TEST_CASE("KCT.util.ChebyshevPolynomialsExplicit.valueAt",
          "[polynomial][chebyshevevaluate][NOPRINT][NOVIZ]")
{
    // First define Chebyshev Polynomial of degree 110
    util::ChebyshevPolynomialsExplicit cp8(
        8, 10, 110, true, 1); // Polynomial of degree 8 that will report eight values
    //    for(int i = 0; i != 20; i++)
    //    {
    //        cp110.printChebyshevPolynomial(i);
    //    }

    double* values = new double[8];
    // Check with https://www.wolframalpha.com/input/?i=+chebyshevT%28110%2Cx%29
    cp8.valuesAt(11.0, values);
    // First we have to evaluate what the time t=11 translates on the domain [-1,1] so
    // ((11-10)/50)-1=-0.99
    double x = -0.98;
    REQUIRE(std::abs(values[0] - x) < 0.000001);
    REQUIRE(std::abs(values[1] - (2 * x * x - 1)) < 0.000001);
    REQUIRE(std::abs(values[2] - (4 * x * x * x - 3 * x)) < 0.000001);
    REQUIRE(std::abs(values[3] - (8 * x * x * x * x - 8 * x * x + 1)) < 0.000001);
    REQUIRE(std::abs(values[4] - (16 * std::pow(x, 5.0) - 20 * std::pow(x, 3.0) + 5 * x))
            < 0.000001);
    REQUIRE(std::abs(values[5]
                     - (32 * std::pow(x, 6.0) - 48 * std::pow(x, 4.0) + 18 * std::pow(x, 2.0) - 1))
            < 0.000001);
    REQUIRE(
        std::abs(values[6]
                 - (64 * std::pow(x, 7.0) - 112 * std::pow(x, 5.0) + 56 * std::pow(x, 3.0) - 7 * x))
        < 0.000001);
    REQUIRE(std::abs(values[7]
                     - (128 * std::pow(x, 8.0) - 256 * std::pow(x, 6.0) + 160 * std::pow(x, 4.0)
                        - 32 * std::pow(x, 2.0) + 1))
            < 0.000001);
    delete[] values;
}

TEST_CASE("KCT.util.ChebyshevPolynomialsExplicit.printChebyshevPolynomial",
          "[polynomial][chebyshev][PRINT][NOVIZ]")
{
    // First define Chebyshev Polynomial of degree 110
    util::ChebyshevPolynomialsExplicit cp110(110, -1, 1, 0);
    for(int i = 0; i != 20; i++)
    {
        cp110.printChebyshevPolynomial(i);
    }
}
