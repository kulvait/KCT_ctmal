// Logging
#include "PLOG/PlogSetup.h"

// Standard libs
#include <iostream>

// Internal libs
#include "MATRIX/Matrix.hpp"
#include "MATRIX/RQFactorization.hpp"
#include "catch.hpp"
#include "stringFormatter.h"
/**First test is simple, just computing Sebastian least squares problem from the excercises.
 *
 */
using namespace CTL;
using namespace CTL::matrix;

/**Test if the procedure will find the correct decomposition
 *
 *A =
 *
 *     1     2     3     5
 *     5    12     1     0
 *     1    -1     1    -3
 *
 *R =
 *   -4.5851   -1.9732    3.7528
 *         0  -12.9228    1.7321
 *         0         0   -3.4641
 *
 *Q =
 *   -0.2712    0.1831   -0.8406   -0.4316
 *   -0.4256   -0.8899   -0.1161    0.1161
 *   -0.2887    0.2887   -0.2887    0.8660
*
*Tested in MATLAB
A = [1 2 3 5; 5 12 1 0; 1 -1 1 -3]

[m n]=size(A);
if m>n
    error('RQ: Number of rows must be smaller than column');
end

[Q R]=qr(flipud(A).');
R=flipud(R.');
R(:,1:m)=R(:,m:-1:1);
Q=Q.';
Q(1:m,:)=Q(m:-1:1,:);
[R, Q] = rq(A);
*
 */

TEST_CASE("TEST: Test RQ factorization on particular example", "RQ factorization")
{

    // Testing if the matrix norm works well
    matrix::Matrix A_(3, 4, { 1, 2, 3, 5, 5, 12, 1, 0, 1, -1, 1, -3 });
    matrix::Matrix R_(3,3,{ -4.5851, -1.9732, 3.7528, 0, -12.9228, 1.7321, 0, 0, -3.4641 });
    matrix::Matrix Q_(3,4,{ -0.2712, 0.1831, -0.8406, -0.4316, -0.4256, -0.8899, -0.1161, 0.1161,
                              -0.2887, 0.2887, -0.2887, 0.8660 });
    std::shared_ptr<matrix::Matrix> A = std::make_shared<matrix::Matrix>(A_);
    std::shared_ptr<matrix::Matrix> R = std::make_shared<matrix::Matrix>(-1.0 * R_);
    std::shared_ptr<matrix::Matrix> Q = std::make_shared<matrix::Matrix>(-1.0 * Q_);
    matrix::RQFactorization rq;
    rq.factorize(A);
    auto R_computed = rq.getRMatrix();
    auto Q_computed = rq.getQMatrix();
    std::cout << A->toString("A");
    std::cout << "Matrix R computed:" << std::endl;
    std::cout << R_computed->toString("R");
    std::cout << "Matrix R expected:" << std::endl;
    std::cout << R->toString("R");
    std::cout << "Matrix Q computed:" << std::endl;
    std::cout << Q_computed->toString("Q");
    std::cout << "Matrix Q expected:" << std::endl;
    std::cout << Q->toString("Q");
    REQUIRE((*R - *R_computed).norm() < 1e-3);
    REQUIRE((*Q - *Q_computed).norm() < 1e-3);
}
