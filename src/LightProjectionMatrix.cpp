#include "MATRIX/LightProjectionMatrix.hpp"

namespace KCT {
namespace matrix {
    // matrix object will be square with the third row unit vector towards detector
    LightProjectionMatrix::LightProjectionMatrix(const double* sourcePosition,
                                                 const double* vector9)
    {
        std::copy(sourcePosition, sourcePosition + 3, std::begin(this->source));
        std::copy(vector9, vector9 + 9, std::begin(this->matrix));
        std::array<double, 3> n;
        std::array<double, 3> v0 = multiplyVectorByConstant(
            source, -1.0); // Coordinates of word point (0.0,0.0,0.0) shifted by source
        std::copy(std::begin(this->matrix) + 6, std::begin(this->matrix) + 9, std::begin(n));
        double matrixMultiplicationFactor = 1.0 / vectorNorm<3>(n);
        if(vectorDotProduct<3>(n, v0) < 0)
        {
            matrixMultiplicationFactor *= -1.0;
        }
        this->matrix = multiplyVectorByConstant<9>(this->matrix, matrixMultiplicationFactor);
    }

    LightProjectionMatrix::LightProjectionMatrix(const ProjectionMatrix& pm)
        : LightProjectionMatrix(std::begin(pm.sourcePosition()),
                                std::begin(pm.projectionMatrixAsVector9()))
    {
    }

    std::array<double, 3> LightProjectionMatrix::sourcePosition() const { return source; }
    double LightProjectionMatrix::pixelSkew() const
    {
        std::array<double, 3> VX = directionVectorVX();
        std::array<double, 3> VY = directionVectorVY();
        VX = normalizeVector<3>(VX);
        VY = normalizeVector<3>(VY);
        return vectorDotProduct<3>(VX, VY);
    }

    std::array<double, 3> LightProjectionMatrix::directionVectorVN() const
    {
        std::array<double, 3> v;
        directionVectorVN(std::begin(v));
        return v;
    }

    std::array<double, 3> LightProjectionMatrix::directionVectorVX() const
    {
        std::array<double, 3> v;
        directionVectorVX(std::begin(v));
        return v;
    }

    std::array<double, 3> LightProjectionMatrix::directionVectorVY() const
    {
        std::array<double, 3> v;
        directionVectorVY(std::begin(v));
        return v;
    }
    std::array<double, 2> LightProjectionMatrix::principalRayProjection() const
    {
        std::array<double, 2> v;
        principalRayProjection(std::begin(v));
        return v;
    }
    std::array<double, 3> LightProjectionMatrix::normalToDetector() const
    {
        std::array<double, 3> v;
        normalToDetector(std::begin(v));
        return v;
    }

    std::array<double, 3> LightProjectionMatrix::directionToPosition(const double pi,
                                                                     const double pj) const
    {
        std::array<double, 3> v;
        directionToPosition(pi, pj, std::begin(v));
        return v;
    }

    void LightProjectionMatrix::project(
        const double x, const double y, const double z, double* pi, double* pj) const
    {
        project0(x - source[0], y - source[1], z - source[2], pi, pj);
    }

    void LightProjectionMatrix::project(
        const float x, const float y, const float z, float* pi, float* pj) const
    {
        project0(x - source[0], y - source[1], z - source[2], pi, pj);
    }
    void LightProjectionMatrix::project0(
        const double X0, const double X1, const double X2, double* pi, double* pj) const
    {
        std::array<double, 3> X{ X0, X1, X2 };
        std::array<double, 3> VX = directionVectorVX();
        std::array<double, 3> VY = directionVectorVY();
        std::array<double, 3> VN = directionVectorVN();
        std::array<double, 2> P0 = principalRayProjection();
        double PX = vectorDotProduct<3>(X, VX);
        double PY = vectorDotProduct<3>(X, VY);
        double PN = vectorDotProduct<3>(X, VN);
        (*pi) = P0[0] + PX / PN;
        (*pj) = P0[1] + PY / PN;
    }

    void LightProjectionMatrix::project0(
        const float X0, const float X1, const float X2, float* pi, float* pj) const
    {
        double pid, pjd;
        project0(X0, X1, X2, &pid, &pjd);
        (*pi) = pid;
        (*pj) = pjd;
    }

    std::array<double, 2> LightProjectionMatrix::focalLength() const
    {
        std::array<double, 2> v;
        focalLength(std::begin(v));
        return v;
    }
    std::array<double, 2> LightProjectionMatrix::pixelSizes(const double sourceToDetector) const
    {
        std::array<double, 2> v;
        pixelSizes(sourceToDetector, std::begin(v));
        return v;
    }
    double LightProjectionMatrix::sourceToDetectorFromPX(const double PX) const
    {
        // PX is pixel size in X direction
        std::array<double, 3> VX = directionVectorVX();
        double fx = vectorNorm<3>(VX);
        return fx * PX;
    }

    double LightProjectionMatrix::sourceToDetectorFromPY(const double PY) const
    {
        // PY is pixel size in Y direction
        std::array<double, 3> VY = directionVectorVY();
        double fy = vectorNorm<3>(VY);
        return fy * PY;
    }

    void LightProjectionMatrix::sourcePosition(double* vector3) const
    {
        std::copy(std::begin(this->source), std::begin(this->source) + 3, vector3);
    }

    void LightProjectionMatrix::directionVectorVN(double* vector3) const
    {
        std::copy(std::begin(this->matrix) + 6, std::begin(this->matrix) + 9, vector3);
    }

    void LightProjectionMatrix::directionVectorVX(double* vector3) const
    {
        std::array<double, 3> VN = directionVectorVN();
        std::array<double, 3> VX;
        std::copy(std::begin(this->matrix), std::begin(this->matrix) + 3, std::begin(VX));
        VX = orthogonalPartOfVectorWithRespectToSecondVector<3>(VX, VN);
        std::copy(std::begin(VX), std::begin(VX) + 3, vector3);
    }

    void LightProjectionMatrix::directionVectorVY(double* vector3) const
    {
        std::array<double, 3> VN = directionVectorVN();
        std::array<double, 3> VY;
        std::copy(std::begin(this->matrix) + 3, std::begin(this->matrix) + 6, std::begin(VY));
        VY = orthogonalPartOfVectorWithRespectToSecondVector<3>(VY, VN);
        std::copy(std::begin(VY), std::begin(VY) + 3, vector3);
    }

    void LightProjectionMatrix::normalToDetector(double* vector3) const
    {
        std::array<double, 3> VN = directionVectorVN();
        VN = multiplyVectorByConstant<3>(VN, -1.0);
        std::copy(std::begin(VN), std::begin(VN) + 3, vector3);
    }

    void LightProjectionMatrix::principalRayProjection(double* vector2) const
    {
        std::array<double, 3> VN = directionVectorVN();
        std::array<double, 3> PX, PY;
        std::copy(std::begin(matrix), std::begin(matrix) + 3, std::begin(PX));
        std::copy(std::begin(matrix) + 3, std::begin(matrix) + 6, std::begin(PY));
        vector2[0] = vectorDotProduct<3>(PX, VN);
        vector2[1] = vectorDotProduct<3>(PY, VN);
    }

    void LightProjectionMatrix::focalLength(double* vector2) const
    {
        std::array<double, 3> VX = directionVectorVX();
        std::array<double, 3> VY = directionVectorVY();
        vector2[0] = vectorNorm<3>(VX);
        vector2[1] = vectorNorm<3>(VY);
    }

    void LightProjectionMatrix::pixelSizes(const double sourceToDetector, double* vector2) const
    {

        std::array<double, 3> VX = directionVectorVX();
        std::array<double, 3> VY = directionVectorVY();
        vector2[0] = sourceToDetector / vectorNorm<3>(VX);
        vector2[1] = sourceToDetector / vectorNorm<3>(VY);
    }
    void LightProjectionMatrix::directionToPosition(const double pi,
                                                    const double pj,
                                                    double* vector3) const
    {
        std::array<double, 3> VX = directionVectorVX();
        double fx2 = vectorDotProduct<3>(VX, VX);
        std::array<double, 3> VY = directionVectorVY();
        double fy2 = vectorDotProduct<3>(VY, VY);
        std::array<double, 3> VN = directionVectorVN();
        std::array<double, 2> P0 = principalRayProjection();
        VX = multiplyVectorByConstant<3>(VX, (pi - P0[0]) / fx2);
        VY = multiplyVectorByConstant<3>(VY, (pj - P0[1]) / fy2);
        VN = vectorSum<3>(VN, VX);
        VN = vectorSum<3>(VN, VY);
        std::copy(std::begin(VN), std::begin(VN) + 3, vector3);
    }

    void LightProjectionMatrix::projectionMatrixAsVector12(double* vector12) const
    {
        std::array<double, 3> P1, P2, P3;

        std::copy(std::begin(matrix), std::begin(matrix) + 3, std::begin(P1));
        std::copy(std::begin(matrix) + 3, std::begin(matrix) + 6, std::begin(P2));
        std::copy(std::begin(matrix) + 6, std::begin(matrix) + 9, std::begin(P3));

        std::copy(std::begin(matrix), std::begin(matrix) + 3, vector12);
        vector12[3] = -vectorDotProduct<3>(P1, source);
        std::copy(std::begin(matrix) + 3, std::begin(matrix) + 6, vector12 + 4);
        vector12[7] = -vectorDotProduct<3>(P2, source);
        std::copy(std::begin(matrix) + 6, std::begin(matrix) + 9, vector12 + 8);
        vector12[11] = -vectorDotProduct<3>(P3, source);
    }

    void LightProjectionMatrix::projectionMatrixAsVector9(double* vector9) const
    {
        std::copy(std::begin(matrix), std::begin(matrix) + 9, vector9);
    }

    void LightProjectionMatrix::inverseProjectionMatrixAsVector16(double* vector16) const
    {
        std::array<double, 16> extendedCameraMatrix;
        projectionMatrixAsVector12(std::begin(extendedCameraMatrix));
        extendedCameraMatrix = multiplyVectorByConstant<16>(
            extendedCameraMatrix,
            1 / std::abs(extendedCameraMatrix[11])); // Scaling to have same representation
        std::copy(std::begin(source), std::begin(source) + 3,
                  std::begin(extendedCameraMatrix) + 12);
        extendedCameraMatrix[15] = 1.0;
        // Inverse computation
        lapack_int n = 4;
        lapack_int inf;
        std::array<lapack_int, 4> ipiv; // Array, size at least max(1,min(m, n)). Contains the pivot
                                        // indices; for 1 ≤i≤ min(m, n), row i was interchanged with
                                        // row ipiv(i).
        inf = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, n, n, std::begin(extendedCameraMatrix), n,
                             std::begin(ipiv)); // LU factorization of matrix
        if(inf != 0)
        {
            std::ostringstream matrixStream;
            matrixStream << "Extended Camera Matrix (4x4):\n";

            for(int i = 0; i < 4; ++i)
            {
                for(int j = 0; j < 4; ++j)
                {
                    matrixStream << std::setw(12) << extendedCameraMatrix[i * 4 + j] << " ";
                }
                matrixStream << "\n";
            }

            std::string MATRIX = matrixStream.str();
            LOGE << io::xprintf("Problem handling maxtrix %s", MATRIX.c_str());
            std::string ERR;
            if(inf < 0)
            {
                ERR = io::xprintf("Parameter %d had ilegal value.", -inf);
            } else
            {
                ERR = io::xprintf(
                    "The factorization has been completed, but U is exactly singular, U_%d,%d=0.",
                    inf, inf);
            }
            KCTERR(ERR);
        }
        inf = LAPACKE_dgetri(LAPACK_ROW_MAJOR, n, std::begin(extendedCameraMatrix), n,
                             std::begin(ipiv));
        if(inf != 0)
        {
            std::string ERR;
            if(inf < 0)
            {
                ERR = io::xprintf("Parameter %d had ilegal value.", -inf);
            } else
            {
                ERR = io::xprintf(
                    "The factorization has been completed, but U is exactly singular, U_%d,%d=0.",
                    inf, inf);
            }
            LOGE << ERR;
            throw std::runtime_error(ERR);
        }
        std::copy(std::begin(extendedCameraMatrix), std::end(extendedCameraMatrix), vector16);
    }

    void LightProjectionMatrix::inverseProjectionMatrixAsVector9(double* vector9) const
    {
        /*
Implementation for zero pixelSkew
std::array<double, 3> VX = directionVectorVX();
std::array<double, 3> VY = directionVectorVY();
std::array<double, 3> VXn = normalizeVector<3>(VX);
std::array<double, 3> VYn = normalizeVector<3>(VY);
double pi = vectorNorm<3>(VX);
double pj = vectorNorm<3>(VY);
std::array<double, 3> VN = directionVectorVN();
VX = multiplyVectorByConstant<3>(VXn, 1 / pi);
VY = multiplyVectorByConstant<3>(VYn, 1 / pj);
VXn = multiplyVectorByConstant<3>(VXn, -1.0);
VYn = multiplyVectorByConstant<3>(VYn, -1.0);
VN = vectorSum<3>(VN, VXn);
VN = vectorSum<3>(VN, VYn);
vector9[0] = VX[0];
vector9[1] = VY[0];
vector9[2] = VN[0];
vector9[3] = VX[1];
vector9[4] = VY[1];
vector9[5] = VN[1];
vector9[6] = VX[2];
vector9[7] = VY[2];
vector9[8] = VN[2];
        */
        std::array<double, 9> invertedCameraMatrix = matrix;
        lapack_int n = 3;
        lapack_int inf;
        std::array<lapack_int, 3> ipiv; // Array, size at least max(1,min(m, n)). Contains the pivot
                                        // indices; for 1 ≤i≤ min(m, n), row i was interchanged with
                                        // row ipiv(i).
        inf = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, n, n, std::begin(invertedCameraMatrix), n,
                             std::begin(ipiv)); // LU factorization of matrix
        if(inf != 0)
        {
            std::string ERR;
            if(inf < 0)
            {
                ERR = io::xprintf("Parameter %d had ilegal value.", -inf);
            } else
            {
                ERR = io::xprintf(
                    "The factorization has been completed, but U is exactly singular, U_%d,%d=0.",
                    inf, inf);
            }
            LOGE << ERR;
            throw std::runtime_error(ERR);
        }
        inf = LAPACKE_dgetri(LAPACK_ROW_MAJOR, n, std::begin(invertedCameraMatrix), n,
                             std::begin(ipiv));
        if(inf != 0)
        {
            std::string ERR;
            if(inf < 0)
            {
                ERR = io::xprintf("Parameter %d had ilegal value.", -inf);
            } else
            {
                ERR = io::xprintf(
                    "The factorization has been completed, but U is exactly singular, U_%d,%d=0.",
                    inf, inf);
            }
            LOGE << ERR;
            throw std::runtime_error(ERR);
        }
        std::copy(std::begin(invertedCameraMatrix), std::end(invertedCameraMatrix), vector9);
    }

    std::array<double, 12> LightProjectionMatrix::projectionMatrixAsVector12() const
    {
        std::array<double, 12> v;
        projectionMatrixAsVector12(std::begin(v));
        return v;
    }
    std::array<double, 9> LightProjectionMatrix::projectionMatrixAsVector9() const
    {

        std::array<double, 9> v;
        projectionMatrixAsVector9(std::begin(v));
        return v;
    }
    std::array<double, 16> LightProjectionMatrix::inverseProjectionMatrixAsVector16() const
    {

        std::array<double, 16> v;
        inverseProjectionMatrixAsVector16(std::begin(v));
        return v;
    }
    std::array<double, 9> LightProjectionMatrix::inverseProjectionMatrixAsVector9() const
    {
        std::array<double, 9> v;
        inverseProjectionMatrixAsVector9(std::begin(v));
        return v;
    }
} // namespace matrix
} // namespace KCT
