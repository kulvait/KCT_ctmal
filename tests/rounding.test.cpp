#include "catch2/catch.hpp"
#include <plog/Log.h>

// Internal libs
#include "stringFormatter.h"
#include "MATH/round.h"

using namespace CTL;

///Testing correctness of legendre polynomials and its derivatives

TEST_CASE("CTL.math.round.h", "[rounding]")
{
REQUIRE(math::roundLow(float(-3.5))==float(-4));
    REQUIRE(math::roundLow(float(-2.5))==float(-3));
    REQUIRE(math::roundLow(float(-1.5))==float(-2));
    REQUIRE(math::roundLow(float(-0.5))==float(-1));
    REQUIRE(math::roundLow(float(0.5))==float(0));
    REQUIRE(math::roundLow(float(1.5))==float(1));
    REQUIRE(math::roundLow(float(2.5))==float(2));
    REQUIRE(math::roundLow(float(3.5))==float(3));
}

