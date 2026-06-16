// Logging
#include "PLOG/PlogSetup.h"

// Standard libs
#include <cmath>
#include <iostream>

// Internal libs
#include "DEN/DenAsyncFrame2DWritter.hpp"
#include "DEN/DenGeometry2DParallelReader.hpp"
#include "FrameMemoryViewer2D.hpp"
#include "MATRIX/LightProjectionMatrix.hpp"
#include "MATRIX/Matrix.hpp"
#include "MATRIX/ProjectionMatrix.hpp"
#include "MATRIX/utils.hpp"
#include "PROG/Program.hpp"
#include "PROG/RunTimeInfo.hpp"
#include "catch.hpp"
#include "helpers.test.hpp"
#include "stringFormatter.h"

/**First test is simple, just computing Sebastian least squares problem from the excercises.
 *
 */
using namespace KCT;
using namespace KCT::util;
using namespace KCT::matrix;

TEST_CASE("Geometry2DParallelCameraMatrix", "[noprint]")
{
    // Testing if the matrix norm works well
    util::RunTimeInfo rti;
    std::string pth = rti.getExecutableDirectoryPath();
    io::DenGeometry2DParallelReader dpr(
        io::xprintf("%s/../tests/testFiles/redExp_pbct_CM.den", pth.c_str()));
    double pixel_sizex = 0.001064;
    REQUIRE(dpr.count() == 5001);
    for(uint32_t k = 0; k != dpr.count(); k++)
    {
        geometry::Geometry2DParallel geom = dpr.readGeometry(k);
        std::array<double, 3> VX = geom.directionVectorVX();
        double pixelSpacing = geom.pixelSpacing();
        REQUIRE(almostEqual(pixelSpacing, pixel_sizex));
        REQUIRE(almostEqual(vectorNorm(VX), pixel_sizex));
    }
}
