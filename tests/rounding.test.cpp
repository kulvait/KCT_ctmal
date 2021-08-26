#include "catch.hpp"
#include <plog/Log.h>

// Internal libs
#include "MATH/round.h"
#include "stringFormatter.h"

using namespace KCT;

/// Testing correctness of legendre polynomials and its derivatives

TEST_CASE("KCT.math.round.h", "[rounding]")
{
    REQUIRE(math::roundLow(float(-3.5)) == float(-4));
    REQUIRE(math::roundLow(float(-2.5)) == float(-3));
    REQUIRE(math::roundLow(float(-1.5)) == float(-2));
    REQUIRE(math::roundLow(float(-0.5)) == float(-1));
    REQUIRE(math::roundLow(float(0.5)) == float(0));
    REQUIRE(math::roundLow(float(1.5)) == float(1));
    REQUIRE(math::roundLow(float(2.5)) == float(2));
    REQUIRE(math::roundLow(float(3.5)) == float(3));

    REQUIRE(math::lroundLow(float(-3.5)) == long(-4));
    REQUIRE(math::lroundLow(float(-2.5)) == long(-3));
    REQUIRE(math::lroundLow(float(-1.5)) == long(-2));
    REQUIRE(math::lroundLow(float(-0.5)) == long(-1));
    REQUIRE(math::lroundLow(float(0.5)) == long(0));
    REQUIRE(math::lroundLow(float(1.5)) == long(1));
    REQUIRE(math::lroundLow(float(2.5)) == long(2));
    REQUIRE(math::lroundLow(float(3.5)) == long(3));
}
