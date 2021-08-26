// Logging
#include "PLOG/PlogSetup.h"

// Standard libs
#include <cmath>
#include <iostream>

// Internal libs
#include "DEN/DenAsyncFrame2DWritter.hpp"
#include "DEN/DenProjectionMatrixReader.hpp"
#include "FrameMemoryViewer2D.hpp"
#include "MATRIX/LightProjectionMatrix.hpp"
#include "MATRIX/Matrix.hpp"
#include "MATRIX/ProjectionMatrix.hpp"
#include "MATRIX/utils.hpp"
#include "PROG/Program.hpp"
#include "PROG/RunTimeInfo.hpp"
#include "catch.hpp"
#include "stringFormatter.h"

/**First test is simple, just computing Sebastian least squares problem from the excercises.
 *
 */
using namespace KCT;
using namespace KCT::util;
using namespace KCT::matrix;

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
    util::RunTimeInfo rti;
    std::string pth = rti.getExecutableDirectoryPath();
    io::DenProjectionMatrixReader dpr(io::xprintf("%s/../tests/phantom.matrices", pth.c_str()));
    ProjectionMatrix pm = dpr.readMatrix(10);
    ProjectionMatrix pm_shift = pm.shiftDetectorOrigin(2.0, 1.0);
    LOGD << std::endl << "Printing original and shifted matrices:";
    LOGD << std::endl << pm.toString();
    LOGD << std::endl << pm_shift.toString();
}

bool almostEqual(double x, double y, double tol = 1e-10)
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

template <int N>
bool almostEqual(std::array<double, N> x, std::array<double, N> y, double tol = 1e-10)
{
    std::array<double, N> minusx = multiplyVectorByConstant<N>(x, -1.0);
    std::array<double, N> vdf = vectorSum<N>(minusx, y);
    double n = vectorNorm<N>(vdf);
    if(n < tol)
    {
        return true;
    } else
    {
        return false;
    }
}

template <int N>
bool almostEqualRelative(std::array<double, N> x, std::array<double, N> y, double tol = 1e-10)
{
    std::array<double, N> minusx = multiplyVectorByConstant<N>(x, -1.0);
    std::array<double, N> vdf = vectorSum<N>(minusx, y);
    double basenorm = vectorNorm<N>(x);
    double n = vectorNorm<N>(vdf);
    if(n / basenorm < tol)
    {
        return true;
    } else
    {
        return false;
    }
}

template <int N>
std::string printVector(std::string name, std::array<double, N> x)
{
    std::string s = io::xprintf("%s: (", name.c_str());
    for(uint32_t i = 0; i < N - 1; i++)
    {
        s += io::xprintf("%f, ", x[i]);
    }
    s += io::xprintf("%f)", x[N - 1]);
    return s;
}

TEST_CASE("LightProjectionMatrix.cpp.testing", "[PM]")
{
    util::Program prog(0, nullptr);
    util::RunTimeInfo rti;
    std::string pth = rti.getExecutableDirectoryPath();
    io::DenProjectionMatrixReader dpr(io::xprintf("%s/../tests/ArtisQ.matrix", pth.c_str()));
    double pixelSpacingX = 0.616, pixelSpacingY = 0.616;

    uint32_t numAngles = dpr.count();
    for(uint32_t k = 0; k != numAngles; k++)
    {
        ProjectionMatrix pm = dpr.readMatrix(k);
        LightProjectionMatrix lpm(pm);
        std::array<double, 3> va, vb;
        std::array<double, 2> pa, pb;
        std::array<double, 12> pea, peb;
        std::array<double, 9> pma, pmb;
        std::array<double, 16> iea, ieb;

        REQUIRE(almostEqual<3>(pm.sourcePosition(), lpm.sourcePosition()));
        pm.sourcePosition(std::begin(va));
        lpm.sourcePosition(std::begin(vb));
        REQUIRE(almostEqual<3>(va, vb));

        REQUIRE(almostEqual<3>(pm.directionVectorVN(), lpm.directionVectorVN()));
        pm.directionVectorVN(std::begin(va));
        lpm.directionVectorVN(std::begin(vb));
        REQUIRE(almostEqual<3>(va, vb));

        pm.directionVectorVX(std::begin(va));
        lpm.directionVectorVX(std::begin(vb));
        REQUIRE(almostEqual<3>(va, vb));
        REQUIRE(almostEqual<3>(pm.directionVectorVX(), lpm.directionVectorVX()));

        REQUIRE(almostEqual<3>(pm.directionVectorVY(), lpm.directionVectorVY()));
        pm.directionVectorVY(std::begin(va));
        lpm.directionVectorVY(std::begin(vb));
        REQUIRE(almostEqual<3>(va, vb));

        REQUIRE(almostEqual(pm.pixelSkew(), lpm.pixelSkew()));

        REQUIRE(almostEqual<3>(pm.normalToDetector(), lpm.normalToDetector()));
        pm.normalToDetector(std::begin(va));
        lpm.normalToDetector(std::begin(vb));
        REQUIRE(almostEqual<3>(va, vb));

        pm.principalRayProjection(std::begin(pa));
        lpm.principalRayProjection(std::begin(pb));
        REQUIRE(almostEqual<2>(pa, pb));
        REQUIRE(almostEqual<2>(pm.principalRayProjection(), lpm.principalRayProjection()));

        REQUIRE(almostEqual<2>(pm.focalLength(), lpm.focalLength()));
        pm.focalLength(std::begin(pa));
        lpm.focalLength(std::begin(pb));
        REQUIRE(almostEqual<2>(pa, pb));

        double sourceToDetector = 1200;
        REQUIRE(almostEqual<2>(pm.pixelSizes(sourceToDetector), lpm.pixelSizes(sourceToDetector)));
        pm.pixelSizes(sourceToDetector, std::begin(pa));
        lpm.pixelSizes(sourceToDetector, std::begin(pb));
        REQUIRE(almostEqual<2>(pa, pb));
        // std::cout << "Iter " << k << " pixel dimensions " << printVector<2>("PA", pa) <<
        // std::endl;
        REQUIRE(almostEqual(pa[0], 0.616, 0.01));
        REQUIRE(almostEqual(pa[1], 0.616, 0.01));
        REQUIRE(almostEqual(pb[0], 0.616, 0.01));
        REQUIRE(almostEqual(pb[1], 0.616, 0.01));

        REQUIRE(almostEqual(pm.sourceToDetectorFromPX(0.616), lpm.sourceToDetectorFromPX(0.616)));
        REQUIRE(almostEqual(pm.sourceToDetectorFromPX(0.616), 1200, 10.0));

        for(double pi = -1000; pi < 1000; pi += 100)
        {
            for(double pj = -1000; pj < 1000; pj += 100)
            {

                pm.directionToPosition(pi, pj, std::begin(va));
                lpm.directionToPosition(pi, pj, std::begin(vb));
                REQUIRE(almostEqual<3>(va, vb));
                REQUIRE(almostEqual<3>(pm.directionToPosition(pi, pj),
                                       lpm.directionToPosition(pi, pj)));
            }
        }
        for(double x0 = -10; x0 < 10; x0 += 1)
        {
            for(double x1 = -10; x1 < 10; x1 += 1)
            {
                for(double x2 = -10; x2 < 10; x2 += 1)
                {
                    pm.project(x0, x1, x2, &pa[0], &pa[1]);
                    lpm.project(x0, x1, x2, &pb[0], &pb[1]);
                    REQUIRE(almostEqual<2>(pa, pb));
                    if(x0 != 0.0 || x1 != 0.0 || x2 != 0.0)
                    {
                        pm.project0(x0, x1, x2, &pa[0], &pa[1]);
                        lpm.project0(x0, x1, x2, &pb[0], &pb[1]);
                        REQUIRE(almostEqualRelative<2>(pa, pb));
                    }
                }
            }
        }
        pm.projectionMatrixAsVector12(std::begin(pea));
        lpm.projectionMatrixAsVector12(std::begin(peb));
        REQUIRE(almostEqual<12>(normalizeVector<12>(pea), normalizeVector<12>(peb)));
        pm.projectionMatrixAsVector9(std::begin(pma));
        lpm.projectionMatrixAsVector9(std::begin(pmb));
        REQUIRE(almostEqual<9>(normalizeVector<9>(pma), normalizeVector<9>(pmb)));
        pm.inverseProjectionMatrixAsVector9(std::begin(pma));
        lpm.inverseProjectionMatrixAsVector9(std::begin(pmb));
        REQUIRE(almostEqual<9>(normalizeVector<9>(pma), normalizeVector<9>(pmb)));
        pm.inverseProjectionMatrixAsVector16(std::begin(iea));
        lpm.inverseProjectionMatrixAsVector16(std::begin(ieb));
        REQUIRE(almostEqual<16>(normalizeVector<16>(iea), normalizeVector<16>(ieb)));

        double xoveryspacing = pixelSpacingX / pixelSpacingY;
        double yoverxspacing = pixelSpacingY / pixelSpacingX;
        double x1, x2, y1, y2;
        std::array<double, 3> sourcePosition = pm.sourcePosition();
        std::array<double, 3> normalToDetector = pm.normalToDetector();
        std::array<double, 3> tangentToDetector = pm.tangentToDetectorYDirection();
        pm.project(sourcePosition[0] - normalToDetector[0], sourcePosition[1] - normalToDetector[1],
                   sourcePosition[2] - normalToDetector[2], &x1, &y1);
        pm.project(sourcePosition[0] - normalToDetector[0] + tangentToDetector[0],
                   sourcePosition[1] - normalToDetector[1] + tangentToDetector[1],
                   sourcePosition[2] - normalToDetector[2] + tangentToDetector[2], &x2, &y2);
        double scalingFactor
            = (x1 - x2) * (x1 - x2) * xoveryspacing + (y1 - y2) * (y1 - y2) * yoverxspacing;
        va = lpm.directionVectorVX();
        vb = lpm.directionVectorVY();
        double sf = vectorNorm<3>(va) * vectorNorm<3>(vb);
        LOGI << io::xprintf("Scaling factor A %f B %f", scalingFactor, sf);
        REQUIRE(almostEqual(sf, scalingFactor, 1e-3));
    }
}

TEST_CASE("ProjectionMatrix.normalToDetector.synthetic", "Normal to detector")
{
    util::RunTimeInfo rti;
    std::string pth = rti.getExecutableDirectoryPath();
    io::DenProjectionMatrixReader dpr(io::xprintf("%s/../tests/circular.matrix", pth.c_str()));
    double sourceToDetector = 1200;
    double sourceToCenter = 750;
    const double pi = std::acos(-1);
    int numAngles = 248;
    for(int k = 0; k != numAngles; k++)
    {
        std::array<double, 3> sourceLocation
            = { sourceToCenter * std::cos(k * 2 * pi / dpr.count()),
                sourceToCenter * std::sin(k * 2 * pi / dpr.count()), 0.0 };
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
        REQUIRE(std::abs(s_norm - sourceToCenter) < 0.000001);
        LOGD << "Normal to detector by PM" << std::endl
             << io::xprintf("Normal=(%f, %f, %f).\n", n[0], n[1], n[2]);
        std::array<double, 3> n_source = { s[0] / s_norm, s[1] / s_norm, s[2] / s_norm };
        LOGD << "Normal to detector by design as 0 to source normalized" << std::endl
             << io::xprintf("Normal=(%f, %f, %f).\n", n_source[0], n_source[1], n_source[2]);
        std::array<double, 3> n_thirdrow = { -pm.get(2, 0), -pm.get(2, 1), -pm.get(2, 2) };
        n_thirdrow = vecnorm(n_thirdrow);
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
            = pm.directionToPosition(616.0 * 0.5 - 0.5, 480.0 * 0.5 - 0.5);
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
        double sd1 = normdiff<2>(ptx, np) * 0.616;
        double sd2 = normdiff<2>(pty, np) * 0.616;
        LOGD << "Distance of source to detector based on the projections of tangent vectors should "
                "be "
                "1200"
             << std::endl
             << io::xprintf("Distances are (x, y) = (%f, %f).", sd1, sd2) << std::endl;
        REQUIRE(std::abs(sd1 - sourceToDetector) < 0.00001);
        REQUIRE(std::abs(sd2 - sourceToDetector) < 0.00001);
    }
}

/*
===ArtisQ.matrix 248 views, 0.616mmx0.616mm pixels, 616x480 grid===
Camera matrices ArtisQ.matrix are taken from
ArtisQMagdeburg/PM/export/PM20161006.RUN01.matrices
These were used for many of experiments in the past.
It describes calibrated geometry of ArtisQ system for the first sweep of perfusion scan.
The  distance  from  the  source  to  theisocenter  is 749mm and  the  distance  from  source to the
detector is 1198mm. Detector  matrix  consists  of 2464×1920detector cells of the dimensions 0.154
mm×0.154 mm. For  differentprotocols  there  is  applied  tiling  of 2×2 with  merged  pixelsize
merged pixel size 0.308 mm×0.308 mm and 1232x960 pixels. ArtisQ.matrix were generated for 4×4
tiling, pixel size is 0.616 mm×0.616 mm and 616x480 pixels.

===SourceQ4.matrix 248 views, 0.616mmx0.616mm pixels, 616x480 grid===
Camera matrix derived from ArtisQ.matrix.
For each given source position it computes the normal vector n from the x,y projection of the source
position. The source is placed 749mm in this direction. The isocenter is taken to be 0 and z is
taken as axis of the rotation, principial ray goes through 0 and source to detector distance is
1198mm. Matrix is constructed for 4x4 tiling with pixel size 0.616 mm×0.616 mm and 616x480 pixels.
Projected position (0,0) on the detector is on the upper left pixel.

===SourceQ2.matrix 248 views, 0.308mmx0.308mm pixels, 1232x960 grid===
Camera matrix derived from ArtisQ.matrix.
For each given source position it computes the normal vector n from the x,y projection of the source
position. The source is placed 749mm in this direction. The isocenter is taken to be 0 and z is
taken as axis of the rotation, principial ray goes through 0 and source to detector distance is
1198mm. Matrix is constructed for 2x2 tiling with pixel size 0.308 mm×0.308 mm and 1232x960 pixels.
Projected position (0,0) on the detector is on the upper left pixel.

===SourceQ1.matrix 248 views, 0.154mmx0.154mm pixels, 2464x1920 grid===
Camera matrix derived from ArtisQ.matrix.
For each given source position it computes the normal vector n from the x,y projection of the source
position. The source is placed 749mm in this direction. The isocenter is taken to be 0 and z is
taken as axis of the rotation, principial ray goes through 0 and source to detector distance is
1198mm. Matrix is constructed for 2x2 tiling with pixel size 0.154 mm×0.154 mm and 2464x1920 pixels.
Projected position (0,0) on the detector is on the upper left pixel.
*/
TEST_CASE("ProjectionMatrix.normalToDetector.siemens", "Create aligned matrices.")
{
    util::Program prog(0, nullptr);
    util::RunTimeInfo rti;
    std::string pth = rti.getExecutableDirectoryPath();
    io::DenProjectionMatrixReader dpr(io::xprintf("%s/../tests/ArtisQ.matrix", pth.c_str()));
    uint32_t numAngles = dpr.count();
    io::DenAsyncFrame2DWritter<double> Q4(io::xprintf("%s/../tests/SourceQ4.matrix", pth.c_str()),
                                          4, 3, numAngles);
    io::DenAsyncFrame2DWritter<double> Q2(io::xprintf("%s/../tests/SourceQ2.matrix", pth.c_str()),
                                          4, 3, numAngles);
    io::DenAsyncFrame2DWritter<double> Q1(io::xprintf("%s/../tests/SourceQ1.matrix", pth.c_str()),
                                          4, 3, numAngles);
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
    matrix::Matrix A3(3, 3, { 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0 });
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
            = pm.directionToPosition(616.0 * 0.5 - 0.5, 480.0 * 0.5 - 0.5);
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
        alpha = std::atan2(-s[1], -s[0]);
        matrix::Matrix X1(4, 4);
        X1(0, 0) = std::sin(alpha);
        X1(0, 1) = -std::cos(alpha);
        X1(1, 0) = std::cos(alpha);
        X1(1, 1) = std::sin(alpha);
        X1(2, 2) = -1.0;
        X1(3, 3) = 1.0;
        matrix::Matrix PM4 = A3 * A24 * A14 * E * X2 * X1;
        matrix::Matrix PM2 = A3 * A22 * A12 * E * X2 * X1;
        matrix::Matrix PM1 = A3 * A21 * A11 * E * X2 * X1;
        Q4.writeFrame(io::FrameMemoryViewer2D<double>(ProjectionMatrix(PM4).getPtr(), 4, 3), k);
        Q2.writeFrame(io::FrameMemoryViewer2D<double>(ProjectionMatrix(PM2).getPtr(), 4, 3), k);
        Q1.writeFrame(io::FrameMemoryViewer2D<double>(ProjectionMatrix(PM1).getPtr(), 4, 3), k);
    }
    avgSourceIsocenter /= numAngles;
    avgSourceDetector /= numAngles;
    LOGD << "Average source to detector  distance is " << avgSourceDetector;
    LOGD << "Average source to isocenter distance is " << avgSourceIsocenter;
}

/*
===ArtisP.matrix 496 views, 0.308mmx0.308mm pixels, 1232x960 grid===
Camera matrices ArtisP.matrix are taken from Robert Frysch, original name was
matrices_file01_run1.bin
The  distance  from  the source  to  theisocenter  is 749mm and  the  distance  from  source to the
detector is 1198mm. Detector  matrix  consists  of 2464×1920detector cells of the dimensions 0.154
mm×0.154 mm. For differentprotocols  there  is  applied  tiling  of 2×2 with  merged  pixelsize
merged pixel size 0.308 mm×0.308 mm and 1232x960 pixels. ArtisP.matrix were generated for 2×2
tiling, pixel size is 0.308 mm×0.308 mm and 1232x960 pixels.

===SourceP4.matrix 496 views, 0.616mmx0.616mm pixels, 616x480 grid===
Camera matrix derived from ArtisQ.matrix.
For each given source position it computes the normal vector n from the x,y projection of the source
position. The source is placed 749mm in this direction. The isocenter is taken to be 0 and z is
taken as axis of the rotation, principial ray goes through 0 and source to detector distance is
1198mm. Matrix is constructed for 4x4 tiling with pixel size 0.616 mm×0.616 mm and 616x480 pixels.
Projected position (0,0) on the detector is on the upper left pixel.

===SourceP2.matrix 496 views, 0.308mmx0.308mm pixels, 1232x960 grid===
Camera matrix derived from ArtisQ.matrix.
For each given source position it computes the normal vector n from the x,y projection of the source
position. The source is placed 749mm in this direction. The isocenter is taken to be 0 and z is
taken as axis of the rotation, principial ray goes through 0 and source to detector distance is
1198mm. Matrix is constructed for 2x2 tiling with pixel size 0.308 mm×0.308 mm and 1232x960 pixels.
Projected position (0,0) on the detector is on the upper left pixel.

===SourceP1.matrix 496 views, 0.154mmx0.154mm pixels, 2464x1920 grid===
Camera matrix derived from ArtisQ.matrix.
For each given source position it computes the normal vector n from the x,y projection of the source
position. The source is placed 749mm in this direction. The isocenter is taken to be 0 and z is
taken as axis of the rotation, principial ray goes through 0 and source to detector distance is
1198mm. Matrix is constructed for 2x2 tiling with pixel size 0.154 mm×0.154 mm and 2464x1920 pixels.
Projected position (0,0) on the detector is on the upper left pixel.
*/
TEST_CASE("ProjectionMatrix.variousALgorithms.siemens", "Create matrices.")
{
    util::Program prog(0, nullptr);
    util::RunTimeInfo rti;
    std::string pth = rti.getExecutableDirectoryPath();
    io::DenProjectionMatrixReader dpr(io::xprintf("%s/../tests/ArtisP.matrix", pth.c_str()));
    uint32_t numAngles = dpr.count();
    io::DenAsyncFrame2DWritter<double> Q4(io::xprintf("%s/../tests/SourceP4.matrix", pth.c_str()),
                                          4, 3, numAngles);
    io::DenAsyncFrame2DWritter<double> Q2(io::xprintf("%s/../tests/SourceP2.matrix", pth.c_str()),
                                          4, 3, numAngles);
    io::DenAsyncFrame2DWritter<double> Q1(io::xprintf("%s/../tests/SourceP1.matrix", pth.c_str()),
                                          4, 3, numAngles);
    double sourceToDetector = 1198;
    double sourceToCenter = 749;
    double pixel_size_x = 0.308;
    double pixel_size_y = 0.308;
    double detector_dim_x = 1232.0;
    double detector_dim_y = 960.0;
    matrix::Matrix X2 = matrix::Matrix::unitDiagonal(4, 4);
    X2(1, 3) = sourceToCenter;
    matrix::Matrix E(3, 4,
                     { 1, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1 / sourceToDetector, 0.0, 0.0 });
    matrix::Matrix A14(
        3, 3,
        { 1 / (2 * pixel_size_x), 0.0, 0.0, 0.0, 1 / (2 * pixel_size_y), 0.0, 0.0, 0.0, 1.0 });
    matrix::Matrix A24(3, 3,
                       { 1.0, 0.0, (0.5 * detector_dim_x - 1.0) * 0.5, 0.0, 1.0,
                         (0.5 * detector_dim_y - 1.0) * 0.5, 0.0, 0.0, 1.0 });
    matrix::Matrix A12(3, 3,
                       { 1 / pixel_size_x, 0.0, 0.0, 0.0, 1 / pixel_size_y, 0.0, 0.0, 0.0, 1.0 });
    matrix::Matrix A22(3, 3,
                       { 1.0, 0.0, (detector_dim_x - 1.0) * 0.5, 0.0, 1.0,
                         (detector_dim_y - 1.0) * 0.5, 0.0, 0.0, 1.0 });
    matrix::Matrix A11(
        3, 3,
        { 1 / (0.5 * pixel_size_x), 0.0, 0.0, 0.0, 1 / (0.5 * pixel_size_y), 0.0, 0.0, 0.0, 1.0 });
    matrix::Matrix A21(3, 3,
                       { 1.0, 0.0, (2 * detector_dim_x - 1.0) * 0.5, 0.0, 1.0,
                         (2 * detector_dim_y - 1.0) * 0.5, 0.0, 0.0, 1.0 });
    matrix::Matrix A3(3, 3, { 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0 });
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
            = pm.directionToPosition(detector_dim_x * 0.5 - 0.5, detector_dim_y * 0.5 - 0.5);
        n_center = { -n_center[0], -n_center[1], -n_center[2] };
        LOGD << "Normal to detector by projection to the center" << std::endl
             << io::xprintf("Normal=(%f, %f, %f).\n", n_center[0], n_center[1], n_center[2]);
        REQUIRE(normdiff<3>(n, n_center) < 0.025);
        double px, py;
        pm.project(s[0] + n_center[0], s[1] + n_center[1], s[2] + n_center[2], &px, &py);
        REQUIRE(std::abs(px - detector_dim_x * 0.5 + 0.5) < 0.000001);
        REQUIRE(std::abs(py - detector_dim_y * 0.5 + 0.5) < 0.000001);
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
        double sd1 = normdiff<2>(ptx, np) * pixel_size_x;
        double sd2 = normdiff<2>(pty, np) * pixel_size_y;
        avgSourceDetector += sd2;
        LOGD << "Distance of source to detector based on the projections of tangent vectors should "
                "be "
                "1200"
             << std::endl
             << io::xprintf("Distances are (x, y) = (%f, %f).", sd1, sd2) << std::endl;
        REQUIRE(std::abs(sd1 - 1200.0) < 15);
        REQUIRE(std::abs(sd2 - 1200.0) < 15);
        alpha = std::atan2(-s[1], -s[0]);
        matrix::Matrix X1(4, 4);
        X1(0, 0) = std::sin(alpha);
        X1(0, 1) = -std::cos(alpha);
        X1(1, 0) = std::cos(alpha);
        X1(1, 1) = std::sin(alpha);
        X1(2, 2) = -1.0;
        X1(3, 3) = 1.0;
        matrix::Matrix PM4 = A3 * A24 * A14 * E * X2 * X1;
        matrix::Matrix PM2 = A3 * A22 * A12 * E * X2 * X1;
        matrix::Matrix PM1 = A3 * A21 * A11 * E * X2 * X1;
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
    // LOGE << X2.toString("SHIFT");
    matrix::Matrix E(3, 4,
                     { 1, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1 / sourceToDetector, 0.0, 0.0 });
    // LOGE << E.toString("E");
    matrix::Matrix A1(3, 3,
                      { 1 / pixel_size_x, 0.0, 0.0, 0.0, 1 / pixel_size_y, 0.0, 0.0, 0.0, 1.0 });
    matrix::Matrix A2(3, 3,
                      { 1.0, 0.0, (detector_dim_x - 1.0) * 0.5, 0.0, 1.0,
                        (detector_dim_y - 1.0) * 0.5, 0.0, 0.0, 1.0 });
    for(int i = 0; i != numAngles; ++i)
    {
        matrix::Matrix X1(4, 4);
        X1(0, 0) = std::sin(pi + i * 2 * pi / numAngles);
        X1(0, 1) = -std::cos(pi + i * 2 * pi / numAngles);
        X1(1, 0) = std::cos(pi + i * 2 * pi / numAngles);
        X1(1, 1) = std::sin(pi + i * 2 * pi / numAngles);
        X1(2, 2) = -1.0;
        X1(3, 3) = 1.0;
        // LOGE << X1.toString("ROTATION");
        matrix::Matrix PM = A2 * A1 * E * X2 * X1;
        ProjectionMatrix projmat(PM);
        io::FrameMemoryViewer2D<double> f(projmat.getPtr(), 4, 3);
        // std::cout << PM.toString("PM");
        // std::cout << projmat.toString("PM");
        cmw.writeFrame(f, i);
    }
}
