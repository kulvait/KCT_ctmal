// Logging
#include "PLOG/PlogSetup.h"

// Standard libs
#include <iostream>

// Internal libs
#include "DEN/DenAsyncFrame2DWritter.hpp"
#include "DEN/DenProjectionMatrixReader.hpp"
#include "FrameMemoryViewer2D.hpp"
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

TEST_CASE("TEST: Create projection matrices", "Circular trajectory matrices")
{
    // Testing if the matrix norm works well
	int numAngles = 248;//360
    io::DenAsyncFrame2DWritter<double> cmw("/tmp/circular.matrix", 4, 3, numAngles);
    double sourceToDetector = 1200;
    double pixel_size_x = 0.412109375;
    double pixel_size_y = 0.412109375;
    double detector_dim_x = 616.0;
    double detector_dim_y = 480.0;
    double fx = sourceToDetector / pixel_size_x;
    double fy = sourceToDetector / pixel_size_y;
    const double pi = std::acos(-1);
    matrix::Matrix pm(3, 4, { fx, 0.0, 0.0, 0.0, 0.0, fy, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0 });
    LOGE << pm.toString("pm");
    matrix::Matrix shift = matrix::Matrix::unitDiagonal(4, 4);
    shift(2, 3) = sourceToDetector / 2.0;
    LOGE << shift.toString("SHIFT");
    matrix::Matrix detectorShift(
        3, 3, { 1.0, 0.0, detector_dim_x / 2.0, 0.0, 1.0, detector_dim_y / 2.0, 0.0, 0.0, 1.0 });
    for(int i = 0; i != numAngles; ++i)
    {
        matrix::Matrix rotation(4, 4);
        rotation(0, 0) = std::cos(i * 2*pi / numAngles);
        rotation(0, 1) = std::sin(i * 2*pi / numAngles);
        rotation(1, 2) = 1.0;
        rotation(2, 0) = std::sin(i * 2*pi / numAngles);
        rotation(2, 1) = -std::cos(i * 2*pi / numAngles);
        rotation(3, 3) = 1.0;
        LOGE << rotation.toString("ROTATION");
        matrix::Matrix PM = detectorShift * pm * shift * rotation;
        ProjectionMatrix projmat(PM);
        io::FrameMemoryViewer2D<double> f(projmat.getPtr(), 4, 3);
        std::cout << PM.toString("PM");
        std::cout << projmat.toString("PM");
        cmw.writeFrame(f, i);
    }
}
