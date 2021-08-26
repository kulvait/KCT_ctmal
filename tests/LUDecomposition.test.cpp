// Logging
#include "PLOG/PlogSetup.h"

// Standard libs
#include <iostream>

// Internal libs
#include "DEN/DenFileInfo.hpp"
#include "MATRIX/LUDoolittleForm.hpp"
#include "MATRIX/Matrix.hpp"
#include "MATRIX/SquareMatrix.hpp"
#include "catch.hpp"
#include "stringFormatter.h"
/**First test is simple, just computing Sebastian least squares problem from the excercises.
 *
 */
using namespace KCT;
using namespace KCT::util;
using namespace KCT::matrix;

TEST_CASE("TEST: Matrix norm", "Matrix norm properties")
{
    // Testing if the matrix norm works well
    SquareMatrix X(3, { 1, 1, 1, 1, 1, 1, 1, 1, 1 });
    REQUIRE(X.norm() == 3.0);
}

TEST_CASE("TEST: SquareMatrix class, test basic operations.", "SquareMatrix.test.cpp")
{
    // Now testing LU decomposition with partitial pivoting
    Matrix A(8, 4, { 1, 1, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 1, 0,
                     1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1 });
    Matrix b(8, 1, { 20, 20, 10, 10, 20, 20, 10, 10 });
    Matrix p(8, 1, { 18, 18, 8, 8, 18, 18, 8, 8 });
    //LOGD << std::endl << A.toString("A");
    SquareMatrix AA(A.T() * A);
    //LOGD << AA.toString("\nA^TA");
    Matrix AA_result(4, 4, { 3, 1, 1, 0, 1, 3, 0, 1, 1, 0, 3, 1, 0, 1, 1, 3 });
    REQUIRE((AA - AA_result).norm()
            == 0.0); // Polymorphism ... derived class might act as superclass
    Matrix Ab = A.T() * b;
    REQUIRE(Ab(1, 0) == 50);
}

TEST_CASE("TEST: LU decomposition.", "SquareMatrix.test.cpp")
{
    double minimalPivotSize = 0.0001;
    Matrix AA(4, 4, { 3, 1, 1, 0, 1, 3, 0, 1, 1, 0, 3, 1, 0, 1, 1, 3 });
    LUDoolittleForm lu = LUDoolittleForm::LUDecomposeDoolittle(AA, minimalPivotSize);
    SquareMatrix L_result(4,
                          { 1, 0, 0, 0, 0.333333333333333, 1, 0, 0, 0.333333333333333, -0.125, 1, 0,
                            0, 0.375, 0.428571428571429, 1 });
    SquareMatrix U_result(4,
                          { 3, 1, 1, 0, 0, 2.66666666666667, -0.333333333333333, 1, 0, 0, 2.625,
                            1.125, 0, 0, 0, 2.14285714285714 });
    REQUIRE((lu.getLMatrix() - L_result).norm() < 1e-10);
    REQUIRE((lu.getUMatrix() - U_result).norm() < 1e-10);
    REQUIRE(lu.getDeterminant() == 45);
    SECTION("Computing inverse")
    {
        SquareMatrix INV_result(4, { 0.466666666666667, -0.2, -0.2, 0.133333333333333, -0.2,
                                     0.466666666666667, 0.133333333333333, -0.2, -0.2,
                                     0.133333333333333, 0.466666666666667, -0.2, 0.133333333333333,
                                     -0.2, -0.2, 0.466666666666667 });
        // LOGD << "Inverse : " << std::endl << INV.info();
        REQUIRE((lu.inverseMatrix() - INV_result).norm() < 1e-10);
    }
    // LOGE << "Tolerance matrix: " << (L-L_result).info() << ".";

    // Now try to permute the matrix A swapping order of rows. Pivoting should work.
    //	A = Matrix<8,4>({1,1,0,0,0,0,1,1,0,1,0,0,0,0,1,0,1,0,1,0,0,1,0,1,1,0,0,0,0,0,0,1});

    //	LOGD << "L:" << std::endl << L.info();
    //	LOGD << "L_result:" << std::endl << L_result.info();

    //	LOGD << "U:" << std::endl << lu.getUMatrix().info();
    //	LOGD << "PA:" << std::endl << lu.getPAMatrix().info();
    //	std::array<int, 4> P = lu.getPermutation();
    SECTION("Partial pivoting test that row permutaion produces the same decomposition")
    {
        Matrix AA_rowPermutation = SquareMatrix(4, { 0, 1, 1, 3, 1, 0, 3, 1, 1, 3, 0, 1, 3, 1, 1, 0 });
        lu = LUDoolittleForm::LUDecomposeDoolittle(AA_rowPermutation, minimalPivotSize);
        SquareMatrix INV = lu.inverseMatrix();
        SquareMatrix INV_result(4,
            { 0.133333333333333, -0.2, -0.2, 0.466666666666667, -0.2, 0.133333333333333,
              0.466666666666667, -0.2, -0.2, 0.466666666666667, 0.133333333333333, -0.2,
              0.466666666666667, -0.2, -0.2, 0.133333333333333 });
        REQUIRE((INV - INV_result).norm() < 1e-10);
        REQUIRE((lu.getLMatrix() - L_result).norm() < 1e-10);
        REQUIRE((lu.getUMatrix() - U_result).norm() < 1e-10);
        REQUIRE((lu.getPAMatrix() - AA).norm() < 1e-10);
        REQUIRE(lu.getOddSwapParity() == false);
        // Now swap only rows 2 and 3
        Matrix AA_swap = SquareMatrix(4, { 3, 1, 1, 0, 1, 0, 3, 1, 1, 3, 0, 1, 0, 1, 1, 3 });
        lu = LUDoolittleForm::LUDecomposeDoolittle(AA_swap, minimalPivotSize);
        INV = lu.inverseMatrix();
        REQUIRE((Matrix::unitDiagonal(4,4) - AA_swap * INV).norm() < 1e-10);
        REQUIRE((lu.getLog10Determinant() + 1.653213) < 0.00001);
        REQUIRE(lu.getDeterminant() == -45);
    }
}
