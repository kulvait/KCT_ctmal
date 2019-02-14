// Logging
#include "PLOG/PlogSetup.h"

// Standard libs
#include <iostream>

// Internal libs
#include "DEN/DenProjectionMatrixReader.hpp"
#include "MATRIX/Matrix.hpp"
#include "MATRIX/ProjectionMatrix.hpp"
#include "catch.hpp"
#include "stringFormatter.h"
/**First test is simple, just computing Sebastian least squares problem from the excercises.
 *
 */
using namespace CTL;
using namespace CTL::util;
using namespace CTL::matrix;

TEST_CASE("TEST: Print projection matrix", "Projection matrix printing")
{
    // Testing if the matrix norm works well
    io::DenProjectionMatrixReader dpr("../tests/phantom.matrices");
    ProjectionMatrix pm = dpr.readMatrix(10);
    ProjectionMatrix pm_shift = pm.shiftDetectorOrigin(2.0, 1.0);
    LOGE << "Printing original and shifted matrices:";
    std::cout << pm.toString();
    std::cout << pm_shift.toString();
}
