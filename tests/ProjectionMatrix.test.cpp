// Logging
#include "PLOG/PlogSetup.h"

// Standard libs
#include <cmath>
#include <iostream>

// Internal libs
#include "DEN/DenAsyncFrame2DWritter.hpp"
#include "DEN/DenProjectionMatrixReader.hpp"
#include "FrameMemoryViewer2D.hpp"
#include "MATRIX/Matrix.hpp"
#include "MATRIX/ProjectionMatrix.hpp"
#include "PROG/Program.hpp"
#include "PROG/RunTimeInfo.hpp"
#include "catch.hpp"
#include "stringFormatter.h"

/**First test is simple, just computing Sebastian least squares problem from the excercises.
 *
 */
using namespace CTL;
using namespace CTL::util;
using namespace CTL::matrix;

// Helper functions before tests
std::array<double, 3> vecnorm(std::array<double, 3> v)
{
    std::array<double, 3> n;
    double na = std::sqrt(v[1] * v[1] + v[2] * v[2] + v[0] * v[0]);
    n[0] = v[0] / na;
    n[1] = v[1] / na;
    n[2] = v[2] / na;
    return n;
}

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

TEST_CASE("ProjectionMatrix.toString", "[print]")
{
    // Testing if the matrix norm works well
    io::DenProjectionMatrixReader dpr("../tests/phantom.matrices");
    ProjectionMatrix pm = dpr.readMatrix(10);
    ProjectionMatrix pm_shift = pm.shiftDetectorOrigin(2.0, 1.0);
    LOGD << std::endl << "Printing original and shifted matrices:";
    LOGD << std::endl << pm.toString();
    LOGD << std::endl << pm_shift.toString();
}

TEST_CASE("ProjectionMatrix.normalToDetector.synthetic", "Normal to detector")
{
    io::DenProjectionMatrixReader dpr("../tests/circular.matrix");
    const double pi = std::acos(-1);
    int numAngles = 248;
    for(int k = 0; k != numAngles; k++)
    {
        std::array<double, 3> sourceLocation = { -600.0 * std::sin(k * 2 * pi / dpr.count()),
                                                 600.0 * std::cos(k * 2 * pi / dpr.count()), 0.0 };
        LOGD << "Source position computed by design of circular trajectory" << std::endl
             << io::xprintf("Source=(%f, %f, %f).\n", sourceLocation[0], sourceLocation[1],
                            sourceLocation[2]);
        ProjectionMatrix pm = dpr.readMatrix(k);
        std::array<double, 3> n = pm.normalToDetector();
        std::array<double, 3> s = pm.sourcePosition();
        double s_norm = std::sqrt(s[0] * s[0] + s[1] * s[1] + s[2] * s[2]);
        LOGD << "Source position computed by PM" << std::endl
             << io::xprintf("Source=(%f, %f, %f) with norm %f.\n", s[0], s[1], s[2], s_norm);
        REQUIRE(normdiff<3>(s, sourceLocation) < 0.000001);
        REQUIRE(std::abs(s_norm - 600) < 0.000001);
        LOGD << "Normal to detector by PM" << std::endl
             << io::xprintf("Normal=(%f, %f, %f).\n", n[0], n[1], n[2]);
        std::array<double, 3> n_source = { s[0] / s_norm, s[1] / s_norm, s[2] / s_norm };
        LOGD << "Normal to detector by design as 0 to source normalized" << std::endl
             << io::xprintf("Normal=(%f, %f, %f).\n", n_source[0], n_source[1], n_source[2]);
        std::array<double, 3> n_thirdrow = { -pm.get(2, 0), -pm.get(2, 1), -pm.get(2, 2) };
        LOGD << "Normal to detector by third row of the PM" << std::endl
             << io::xprintf("Normal=(%f, %f, %f).\n", n_thirdrow[0], n_thirdrow[1], n_thirdrow[2]);
        REQUIRE(normdiff<3>(n, n_source) < 0.0000001);
        REQUIRE(normdiff<3>(n, n_thirdrow) < 0.0000001);
        double psx, psy;
        pm.project(s[0], s[1], s[2], &psx, &psy);
        LOGD << "Source projection should be nan" << std::endl
             << io::xprintf("Source projection (px,py) = (%f, %f)\n", psx, psy);
        REQUIRE(std::isnan(psx));
        REQUIRE(std::isnan(psy));
        std::array<double, 3> n_center
            = pm.projectedToPosition(616.0 * 0.5 - 0.5, 480.0 * 0.5 - 0.5);
        n_center = { -n_center[0], -n_center[1], -n_center[2] };
        REQUIRE(normdiff<3>(n, n_center) < 0.0000001);
        double px, py;
        pm.project(s[0] + n_center[0], s[1] + n_center[1], s[2] + n_center[2], &px, &py);
        REQUIRE(std::abs(px - 616.0 * 0.5 + 0.5) < 0.000001);
        REQUIRE(std::abs(py - 480.0 * 0.5 + 0.5) < 0.000001);
        std::array<double, 2> np;
        pm.project(s[0] - n[0], s[1] - n[1], s[2] - n[2], &np[0], &np[1]);
        std::array<double, 3> tx = pm.tangentToDetectorXDirection();
        std::array<double, 3> ty = pm.tangentToDetectorYDirection();
        std::array<double, 2> ptx, pty;
        pm.project(s[0] + tx[0], s[1] + tx[1], s[2] + tx[2], &ptx[0], &ptx[1]);
        pm.project(s[0] + ty[0], s[1] + ty[1], s[2] + ty[2], &pty[0], &pty[1]);
        REQUIRE(std::isnan(ptx[0]));
        REQUIRE(std::isnan(ptx[1]));
        REQUIRE(std::isnan(pty[0]));
        REQUIRE(std::isnan(pty[1]));
        LOGD << "Source plus tangents vector projection should be nan" << std::endl
             << io::xprintf("ptx = (%f, %f) pty=(%f, %f)", ptx[0], ptx[1], pty[0], pty[1]);
        pm.project(s[0] - n[0] + tx[0], s[1] - n[1] + tx[1], s[2] - n[2] + tx[2], &ptx[0], &ptx[1]);
        pm.project(s[0] - n[0] + ty[0], s[1] - n[1] + ty[1], s[2] - n[2] + ty[2], &pty[0], &pty[1]);
        double sd1 = normdiff<2>(ptx, np) * 0.412109375;
        double sd2 = normdiff<2>(pty, np) * 0.412109375;
        LOGD << "Distance of source to detector based on the projections of tangent vectors should "
                "be "
                "1200"
             << std::endl
             << io::xprintf("Distances are (x, y) = (%f, %f).", sd1, sd2) << std::endl;
        REQUIRE(std::abs(sd1 - 1200.0) < 0.00001);
        REQUIRE(std::abs(sd2 - 1200.0) < 0.00001);
    }
}

/*
===camera.matrix or ArtisQ.matrix 248 views, 0.616mmx0.616mm pixels, 616x480 grid===
Camera matrices ArtisQ.matrix are taken from
ArtisQMagdeburg/PM/export/PM20161006.RUN01.matrices
These were used for many of experiments in the past.
It describes calibrated geometry of ArtisQ system for the first sweep of perfusion scan.
The  distance  from  the  source  to  theisocenter  is 749mm and  the  distance  from  source to the detector is 1198mm.
Detector  matrix  consists  of 2464×1920detector cells of the dimensions 0.154 mm×0.154 mm. 
For  differentprotocols  there  is  applied  tiling  of 2×2 with  merged  pixelsize merged pixel size 0.308 mm×0.308 mm and 1232x960 pixels.
ArtisQ.matrix were generated for 4×4 tiling, pixel size is 0.616 mm×0.616 mm and 616x480 pixels.

===SourceQ4.matrix 248 views, 0.616mmx0.616mm pixels, 616x480 grid===
Camera matrix derived from ArtisQ.matrix. 
For each given source position it computes the normal vector n from the x,y projection of the source position. The source is placed 749mm in this direction.
The isocenter is taken to be 0 and z is taken as axis of the rotation, principial ray goes through 0 and source to detector distance is 1198mm.
Matrix is constructed for 4x4 tiling with pixel size 0.616 mm×0.616 mm and 616x480 pixels.
Projected position (0,0) on the detector is on the upper left pixel. 

===SourceQ2.matrix 248 views, 0.308mmx0.308mm pixels, 1232x960 grid===
Camera matrix derived from ArtisQ.matrix. 
For each given source position it computes the normal vector n from the x,y projection of the source position. The source is placed 749mm in this direction.
The isocenter is taken to be 0 and z is taken as axis of the rotation, principial ray goes through 0 and source to detector distance is 1198mm.
Matrix is constructed for 2x2 tiling with pixel size 0.308 mm×0.308 mm and 1232x960 pixels.
Projected position (0,0) on the detector is on the upper left pixel. 

===SourceQ1.matrix 248 views, 0.154mmx0.154mm pixels, 2464x1920 grid===
Camera matrix derived from ArtisQ.matrix. 
For each given source position it computes the normal vector n from the x,y projection of the source position. The source is placed 749mm in this direction.
The isocenter is taken to be 0 and z is taken as axis of the rotation, principial ray goes through 0 and source to detector distance is 1198mm.
Matrix is constructed for 2x2 tiling with pixel size 0.154 mm×0.154 mm and 2464x1920 pixels.
Projected position (0,0) on the detector is on the upper left pixel. 
*/
TEST_CASE("ProjectionMatrix.normalToDetector.siemens", "Create aligned matrices.")
{
    util::Program prog(0, nullptr);
    util::RunTimeInfo rti;
    std::string pth = rti.getExecutableDirectoryPath();
    io::DenProjectionMatrixReader dpr(io::xprintf("%s/../tests/camera.matrix", pth.c_str()));
    uint32_t numAngles = dpr.count();
    io::DenAsyncFrame2DWritter<double> Q4(io::xprintf("%s/../tests/SourceQ4.matrix", pth.c_str()), 4, 3, numAngles);
    io::DenAsyncFrame2DWritter<double> Q2(io::xprintf("%s/../tests/SourceQ2.matrix", pth.c_str()), 4, 3, numAngles);
    io::DenAsyncFrame2DWritter<double> Q1(io::xprintf("%s/../tests/SourceQ1.matrix", pth.c_str()), 4, 3, numAngles);
    double sourceToDetector = 1198;
    double sourceToCenter = 749;
    double pixel_size_x = 0.616;
    double pixel_size_y = 0.616;
    double detector_dim_x = 616.0;
    double detector_dim_y = 480.0;
    matrix::Matrix X2 = matrix::Matrix::unitDiagonal(4, 4);
    X2(1, 3) = sourceToCenter;
    matrix::Matrix E(3, 4,
                     { 1, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1 / sourceToDetector, 0.0, 0.0 });
    matrix::Matrix A14(3, 3,
                       { 1 / pixel_size_x, 0.0, 0.0, 0.0, 1 / pixel_size_y, 0.0, 0.0, 0.0, 1.0 });
    matrix::Matrix A24(3, 3,
                       { 1.0, 0.0, (detector_dim_x - 1.0) * 0.5, 0.0, 1.0,
                         (detector_dim_y - 1.0) * 0.5, 0.0, 0.0, 1.0 });
    matrix::Matrix A12(3, 3,
                       { 2 / pixel_size_x, 0.0, 0.0, 0.0, 2 / pixel_size_y, 0.0, 0.0, 0.0, 1.0 });
    matrix::Matrix A22(3, 3,
                       { 1.0, 0.0, (2 * detector_dim_x - 1.0) * 0.5, 0.0, 1.0,
                         (2 * detector_dim_y - 1.0) * 0.5, 0.0, 0.0, 1.0 });
    matrix::Matrix A11(3, 3,
                       { 1 / (0.25 * pixel_size_x), 0.0, 0.0, 0.0, 1 / (0.25 * pixel_size_y), 0.0,
                         0.0, 0.0, 1.0 });
    matrix::Matrix A21(3, 3,
                       { 1.0, 0.0, (4 * detector_dim_x - 1.0) * 0.5, 0.0, 1.0,
                         (4 * detector_dim_y - 1.0) * 0.5, 0.0, 0.0, 1.0 });
    double avgSourceDetector = 0.0;
    double avgSourceIsocenter = 0.0;
    double alpha;
    for(uint32_t k = 0; k != numAngles; k++)
    {
        ProjectionMatrix pm = dpr.readMatrix(k);
        std::array<double, 3> n = pm.normalToDetector();
        std::array<double, 3> s = pm.sourcePosition();
        double s_norm = std::sqrt(s[0] * s[0] + s[1] * s[1] + s[2] * s[2]);
        LOGD << "Source position computed by PM" << std::endl
             << io::xprintf("Source=(%f, %f, %f) with norm %f.\n", s[0], s[1], s[2], s_norm);
        REQUIRE(std::abs(s_norm - 750) < 5);
        LOGD << "Normal to detector by PM" << std::endl
             << io::xprintf("Normal=(%f, %f, %f).\n", n[0], n[1], n[2]);
        std::array<double, 3> n_source = { s[0] / s_norm, s[1] / s_norm, s[2] / s_norm };
        LOGD << "Normal to detector by design as 0 to source normalized" << std::endl
             << io::xprintf("Normal=(%f, %f, %f).\n", n_source[0], n_source[1], n_source[2]);
        std::array<double, 3> n_thirdrow = { -pm.get(2, 0), -pm.get(2, 1), -pm.get(2, 2) };
        n_thirdrow = vecnorm(n_thirdrow);
        LOGD << "Normal to detector by third row of the PM" << std::endl
             << io::xprintf("Normal=(%f, %f, %f).\n", n_thirdrow[0], n_thirdrow[1], n_thirdrow[2]);
        REQUIRE(normdiff<3>(n, n_source) < 0.025);
        REQUIRE(normdiff<3>(n, n_thirdrow) < 0.0000001);
        double psx, psy;
        pm.project(s[0], s[1], s[2], &psx, &psy);
        LOGD << "Source projection should be nan" << std::endl
             << io::xprintf("Source projection (px,py) = (%f, %f)\n", psx, psy);
        avgSourceIsocenter += std::sqrt(s[0] * s[0] + s[1] * s[1]);
        REQUIRE(std::isnan(psx));
        REQUIRE(std::isnan(psy));
        std::array<double, 3> n_center
            = pm.projectedToPosition(616.0 * 0.5 - 0.5, 480.0 * 0.5 - 0.5);
        n_center = { -n_center[0], -n_center[1], -n_center[2] };
        LOGD << "Normal to detector by projection to the center" << std::endl
             << io::xprintf("Normal=(%f, %f, %f).\n", n_center[0], n_center[1], n_center[2]);
        REQUIRE(normdiff<3>(n, n_center) < 0.025);
        double px, py;
        pm.project(s[0] + n_center[0], s[1] + n_center[1], s[2] + n_center[2], &px, &py);
        REQUIRE(std::abs(px - 616.0 * 0.5 + 0.5) < 0.000001);
        REQUIRE(std::abs(py - 480.0 * 0.5 + 0.5) < 0.000001);
        std::array<double, 2> np;
        pm.project(s[0] - n[0], s[1] - n[1], s[2] - n[2], &np[0], &np[1]);
        std::array<double, 3> tx = pm.tangentToDetectorXDirection();
        std::array<double, 3> ty = pm.tangentToDetectorYDirection();
        std::array<double, 2> ptx, pty;
        pm.project(s[0] + tx[0], s[1] + tx[1], s[2] + tx[2], &ptx[0], &ptx[1]);
        pm.project(s[0] + ty[0], s[1] + ty[1], s[2] + ty[2], &pty[0], &pty[1]);
        REQUIRE(std::isnan(ptx[0]));
        REQUIRE(std::isnan(ptx[1]));
        REQUIRE(std::isnan(pty[0]));
        REQUIRE(std::isnan(pty[1]));
        LOGD << "Source plus tangents vector projection should be nan" << std::endl
             << io::xprintf("ptx = (%f, %f) pty=(%f, %f)", ptx[0], ptx[1], pty[0], pty[1]);
        pm.project(s[0] - n[0] + tx[0], s[1] - n[1] + tx[1], s[2] - n[2] + tx[2], &ptx[0], &ptx[1]);
        pm.project(s[0] - n[0] + ty[0], s[1] - n[1] + ty[1], s[2] - n[2] + ty[2], &pty[0], &pty[1]);
        double sd1 = normdiff<2>(ptx, np) * 0.616;
        double sd2 = normdiff<2>(pty, np) * 0.616;
        avgSourceDetector += sd2;
        LOGD << "Distance of source to detector based on the projections of tangent vectors should "
                "be "
                "1200"
             << std::endl
             << io::xprintf("Distances are (x, y) = (%f, %f).", sd1, sd2) << std::endl;
        REQUIRE(std::abs(sd1 - 1200.0) < 15);
        REQUIRE(std::abs(sd2 - 1200.0) < 15);
        alpha = std::atan2(s[1], s[0]);
        matrix::Matrix X1(4, 4);
        X1(0, 0) = -std::sin(alpha);
        X1(0, 1) = std::cos(alpha);
        X1(1, 0) = std::cos(alpha);
        X1(1, 1) = std::sin(alpha);
        X1(2, 2) = -1.0;
        X1(3, 3) = 1.0;
        matrix::Matrix PM4 = A24 * A14 * E * X2 * X1;
        matrix::Matrix PM2 = A22 * A12 * E * X2 * X1;
        matrix::Matrix PM1 = A21 * A11 * E * X2 * X1;
        Q4.writeFrame(io::FrameMemoryViewer2D<double>(ProjectionMatrix(PM4).getPtr(), 4, 3), k);
        Q2.writeFrame(io::FrameMemoryViewer2D<double>(ProjectionMatrix(PM2).getPtr(), 4, 3), k);
        Q1.writeFrame(io::FrameMemoryViewer2D<double>(ProjectionMatrix(PM1).getPtr(), 4, 3), k);
    }
    avgSourceIsocenter /= numAngles;
    avgSourceDetector /= numAngles;
    LOGD << "Average source to detector  distance is " << avgSourceDetector;
    LOGD << "Average source to isocenter distance is " << avgSourceIsocenter;
}

TEST_CASE("CreateTestMatrices", "Circular trajectory matrices")
{
    // Testing if the matrix norm works well
    int numAngles = 248; // 360
    io::DenAsyncFrame2DWritter<double> cmw("/tmp/circular.matrix", 4, 3, numAngles);
    double sourceToDetector = 1200;
    double sourceToCenter = 750;
    double pixel_size_x = 0.616;
    double pixel_size_y = 0.616;
    double detector_dim_x = 616.0;
    double detector_dim_y = 480.0;
    const double pi = std::acos(-1);
    matrix::Matrix X2 = matrix::Matrix::unitDiagonal(4, 4);
    X2(1, 3) = sourceToCenter;
    LOGE << X2.toString("SHIFT");
    matrix::Matrix E(3, 4,
                     { 1, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1 / sourceToDetector, 0.0, 0.0 });
    LOGE << E.toString("E");
    matrix::Matrix A1(3, 3,
                      { 1 / pixel_size_x, 0.0, 0.0, 0.0, 1 / pixel_size_y, 0.0, 0.0, 0.0, 1.0 });
    matrix::Matrix A2(3, 3,
                      { 1.0, 0.0, (detector_dim_x - 1.0) * 0.5, 0.0, 1.0,
                        (detector_dim_y - 1.0) * 0.5, 0.0, 0.0, 1.0 });
    for(int i = 0; i != numAngles; ++i)
    {
        matrix::Matrix X1(4, 4);
        X1(0, 0) = -std::sin(i * 2 * pi / numAngles);
        X1(0, 1) = std::cos(i * 2 * pi / numAngles);
        X1(1, 0) = std::cos(i * 2 * pi / numAngles);
        X1(1, 1) = std::sin(i * 2 * pi / numAngles);
        X1(2, 2) = -1.0;
        X1(3, 3) = 1.0;
        LOGE << X1.toString("ROTATION");
        matrix::Matrix PM = A2 * A1 * E * X2 * X1;
        ProjectionMatrix projmat(PM);
        io::FrameMemoryViewer2D<double> f(projmat.getPtr(), 4, 3);
        std::cout << PM.toString("PM");
        std::cout << projmat.toString("PM");
        cmw.writeFrame(f, i);
    }
}
