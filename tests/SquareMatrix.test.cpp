// Logging
#include "PLOG/PlogSetup.h"

// Standard libs
#include <iostream>

// Internal libs
#include "DEN/DenFileInfo.hpp"
#include "MATRIX/LUDoolittleForm.hpp"
#include "MATRIX/SquareMatrix.hpp"
#include "MATRIX/Matrix.hpp"
#include "catch.hpp"
#include "stringFormatter.h"
/**First test is simple, just computing Sebastian least squares problem from the excercises.
 *
 */
using namespace CTL;
using namespace CTL::util;
using namespace CTL::matrix;

TEST_CASE("TEST: Matrix norm", "Matrix norm properties")
{
    // Testing if the matrix norm works well
    SquareMatrix<3> X({ 1, 1, 1, 1, 1, 1, 1, 1, 1 });
    REQUIRE(X.norm() == 3.0);
}

TEST_CASE("TEST: SquareMatrix class, test basic operations.", "SquareMatrix.test.cpp")
{
    // Now testing LU decomposition with partitial pivoting
    Matrix<8, 4> A({ 1, 1, 0, 0,
			 0, 0, 1, 1,
			 0, 1, 0, 0,
			 0, 0, 1, 0,
                     1, 0, 1, 0,
			 0, 1, 0, 1,
			 1, 0, 0, 0,
			 0, 0, 0, 1 });
    Matrix<8, 1> b({ 20, 20, 10, 10, 20, 20, 10, 10 });
    Matrix<8, 1> p({ 18, 18, 8, 8, 18, 18, 8, 8 });
    // LOGD << "A="<< std::endl << A.info();
    SquareMatrix<4> AA(A.T() * A);
    // LOGD << "A^TA=" << std::endl << AA.info();
    Matrix<4, 4> AA_result({ 3, 1, 1, 0, 1, 3, 0, 1, 1, 0, 3, 1, 0, 1, 1, 3 });
    REQUIRE((AA - AA_result).norm()
            == 0.0); // Polymorphism ... derived class might act as superclass
    Matrix<4, 1> Ab = A.T() * b;
    REQUIRE(Ab(1,0) == 50);
    // LOGD << "A^tp=" << std::endl << Ap.info();
    double minimalPivotSize = 0.0001;
    LUDoolittleForm<4> lu = LUDoolittleForm<4>::LUDecomposeDoolittle(AA, minimalPivotSize);
    SquareMatrix<4> L = lu.getLMatrix();
    SquareMatrix<4> L_result({ 1, 0, 0, 0, 0.333333333333333, 1, 0, 0, 0.333333333333333, -0.125, 1,
                               0, 0, 0.375, 0.428571428571429, 1 });
    SquareMatrix<4> U_result({ 3, 1, 1, 0, 0, 2.66666666666667, -0.333333333333333, 1, 0, 0, 2.625,
                               1.125, 0, 0, 0, 2.14285714285714 });
    REQUIRE((L - L_result).norm() < 1e-10);
    REQUIRE((lu.getUMatrix() - U_result).norm() < 1e-10);
    REQUIRE(lu.getDeterminant() == 45);
    SECTION("Computing inverse")
    {
        auto INV = lu.inverseMatrix();
        SquareMatrix<4> INV_result({ 0.466666666666667, -0.2, -0.2, 0.133333333333333, -0.2,
                                     0.466666666666667, 0.133333333333333, -0.2, -0.2,
                                     0.133333333333333, 0.466666666666667, -0.2, 0.133333333333333,
                                     -0.2, -0.2, 0.466666666666667 });
        // LOGD << "Inverse : " << std::endl << INV.info();
        REQUIRE((INV - INV_result).norm() < 1e-10);
    }
    // LOGE << "Tolerance matrix: " << (L-L_result).info() << ".";

    // Now try to permute the matrix A swapping order of rows. Pivoting should work.
    //	A = Matrix<8,4>({1,1,0,0,0,0,1,1,0,1,0,0,0,0,1,0,1,0,1,0,0,1,0,1,1,0,0,0,0,0,0,1});

    //	LOGD << "L:" << std::endl << L.info();
    //	LOGD << "L_result:" << std::endl << L_result.info();

    //	LOGD << "U:" << std::endl << lu.getUMatrix().info();
    //	LOGD << "PA:" << std::endl << lu.getPAMatrix().info();
    //	std::array<int, 4> P = lu.getPermutation();
    SECTION("Now lets test partial pivoting")
    {
        AA = SquareMatrix<4>({ 0, 1, 1, 3, 1, 0, 3, 1, 1, 3, 0, 1, 3, 1, 1, 0 });
        lu = LUDoolittleForm<4>::LUDecomposeDoolittle(AA, minimalPivotSize);
        SquareMatrix<4> INV = lu.inverseMatrix();
        SquareMatrix<4> INV_result = SquareMatrix<4>(
            { 0.133333333333333, -0.2, -0.2, 0.466666666666667, -0.2, 0.133333333333333,
              0.466666666666667, -0.2, -0.2, 0.466666666666667, 0.133333333333333, -0.2,
              0.466666666666667, -0.2, -0.2, 0.133333333333333 });
        REQUIRE((INV - INV_result).norm() < 1e-10);
        REQUIRE((lu.getLMatrix() - L_result).norm() < 1e-10);
        REQUIRE((lu.getUMatrix() - U_result).norm() < 1e-10);
        REQUIRE((lu.getPAMatrix() - AA_result).norm() < 1e-10);
        REQUIRE(lu.getOddSwapParity() == false);
        // Now swap only rows 2 and 3
        AA = SquareMatrix<4>({ 3, 1, 1, 0, 1, 0, 3, 1, 1, 3, 0, 1, 0, 1, 1, 3 });
        lu = LUDoolittleForm<4>::LUDecomposeDoolittle(AA, minimalPivotSize);
        INV = lu.inverseMatrix();
        REQUIRE((Matrix<4,4>::unitDiagonal() - AA * INV).norm() < 1e-10);
        REQUIRE((lu.getLog10Determinant() + 1.653213) < 0.00001);
        REQUIRE(lu.getDeterminant() == -45);
    }
}
