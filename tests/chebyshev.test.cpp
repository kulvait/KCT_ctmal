#include "catch.hpp"
#include <cmath>
#include <plog/Log.h>

// Internal libs
#include "FUN/ChebyshevPolynomialsExplicit.hpp"
#include "stringFormatter.h"

using namespace CTL;

/// Testing correctness of legendre polynomials and its derivatives

TEST_CASE("CTL.util.ChebyshevPolynomialsExplicit.evaluate",
          "[polynomial][chebyshev][NOPRINT][NOVIZ]")
{
    // First define Chebyshev Polynomial of degree 110
    util::ChebyshevPolynomialsExplicit cp110(110, -1, 1, 0);
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

TEST_CASE("CTL.util.ChebyshevPolynomialsExplicit.printChebyshevPolynomial",
          "[polynomial][chebyshev][PRINT][NOVIZ]")
{
    // First define Chebyshev Polynomial of degree 110
    util::ChebyshevPolynomialsExplicit cp110(110, -1, 1, 0);
    for(int i = 0; i != 20; i++)
    {
        cp110.printChebyshevPolynomial(i);
    }
}
