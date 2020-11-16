#include "MATRIX/ProjectionMatrix.hpp"

namespace CTL {
namespace matrix {

    SquareMatrix ProjectionMatrix::colSubMatrix(int j) const
    {
        if(j < 0 || j > 3)
        {
            std::string errMsg = io::xprintf(
                "Row must be specified in the range [0,3], specified i=%d is out of range.", j);
            LOGE << errMsg;
            throw std::runtime_error(errMsg);
        }
        SquareMatrix out(3);
        for(int i = 0; i != 3; i++)
        {
            for(int jtm = 0; jtm != 4; jtm++)
            {
                if(j > jtm)
                {
                    out(i, jtm) = (*this)(i, jtm);
                } else if(j == jtm)
                {
                    continue;
                } else
                {
                    out(i, jtm - 1) = (*this)(i, jtm);
                }
            }
        }
        return (out);
    }

    double ProjectionMatrix::pixelSkew() const
    {
        std::array<double, 3> VX = directionVectorVX();
        std::array<double, 3> VY = directionVectorVY();
        VX = normalizeVector<3>(VX);
        VY = normalizeVector<3>(VY);
        return vectorDotProduct<3>(VX, VY);
    }

    std::array<double, 3> ProjectionMatrix::directionVectorVN() const
    {
        std::array<double, 3> VN = normalToDetector();
        return matrix::multiplyVectorByConstant<3>(VN, -1.0);
    }

    std::array<double, 3> ProjectionMatrix::directionVectorVX() const
    {
        std::array<double, 3> VN = directionVectorVN();
        std::array<double, 3> VXn = tangentToDetectorXDirection();
        double px, px0, py;
        project0(VN[0] + VXn[0], VN[1] + VXn[1], VN[2] + VXn[2], &px, &py);
        project0(VN[0], VN[1], VN[2], &px0, &py);
        return matrix::multiplyVectorByConstant<3>(VXn, px - px0);
    }
    std::array<double, 3> ProjectionMatrix::directionVectorVY() const
    {

        std::array<double, 3> VN = directionVectorVN();
        std::array<double, 3> VYn = tangentToDetectorYDirection();
        double px, py, py0;
        project0(VN[0] + VYn[0], VN[1] + VYn[1], VN[2] + VYn[2], &px, &py);
        project0(VN[0], VN[1], VN[2], &px, &py0);
        return matrix::multiplyVectorByConstant<3>(VYn, py - py0);
    }

    void ProjectionMatrix::directionVectorVN(double* vector3) const
    {
        std::array<double, 3> V = directionVectorVN();
        std::copy(std::begin(V), std::begin(V) + 3, vector3);
    }

    void ProjectionMatrix::directionVectorVX(double* vector3) const
    {
        std::array<double, 3> V = directionVectorVX();
        std::copy(std::begin(V), std::begin(V) + 3, vector3);
    }
    void ProjectionMatrix::directionVectorVY(double* vector3) const
    {
        std::array<double, 3> V = directionVectorVY();
        std::copy(std::begin(V), std::begin(V) + 3, vector3);
    }

    std::array<double, 3> ProjectionMatrix::tangentToDetectorXDirection() const
    {
        std::array<double, 3> thirdRow = { A[8], A[9], A[10] };
        std::array<double, 3> firstRow = { A[0], A[1], A[2] };
        std::array<double, 3> tangent
            = orthogonalPartOfVectorWithRespectToSecondVector<3>(firstRow, thirdRow);
        return normalizeVector(tangent);
    }

    std::array<double, 3> ProjectionMatrix::tangentToDetectorYDirection() const
    {
        std::array<double, 3> thirdRow = { A[8], A[9], A[10] };
        std::array<double, 3> secondRow = { A[4], A[5], A[6] };
        std::array<double, 3> tangent
            = orthogonalPartOfVectorWithRespectToSecondVector<3>(secondRow, thirdRow);
        return normalizeVector(tangent);
    }

    std::array<double, 3> ProjectionMatrix::sourcePosition() const { return source; }

    /**
     * Normal to the detector. When its added to the source position, its pointing from detector.
     */
    std::array<double, 3> ProjectionMatrix::normalToDetector() const
    {
        std::array<double, 3> normal, tangentxnorm, tangentynorm;
        tangentxnorm = tangentToDetectorXDirection();
        tangentynorm = tangentToDetectorYDirection();
        // Cross product
        normal[0] = tangentxnorm[1] * tangentynorm[2] - tangentxnorm[2] * tangentynorm[1];
        normal[1] = tangentxnorm[2] * tangentynorm[0] - tangentxnorm[0] * tangentynorm[2];
        normal[2] = tangentxnorm[0] * tangentynorm[1] - tangentxnorm[1] * tangentynorm[0];
        // This however will depend on the right or left handedness of the detector
        // Lets see if this is codirectional with the approximate normal from 0 to source position
        if(source[0] * normal[0] + source[1] * normal[1] + source[2] * normal[2] < 0.0)
        {
            normal[0] = -normal[0];
            normal[1] = -normal[1];
            normal[2] = -normal[2];
        }
        // Renormalization
        return (normalizeVector(normal));
    }

    /**
     *Finds the unit vector such that S + (x,y,z,0) is projected to given position on the
     *detector. *For (px,py) center of the detector we should get normal to detector. *Given the
     *structure of the projection matrix being rotation matrix in the first three columns, *the rows
     *should be OG to each other and the third row should describe what projects to (0,0) *before
     *detector shift to upper right corner. Test accuracy of this by calling *normalToDetector(px,
     *py), where px, py are dtetector coordinates.
     *
     * @param px
     * @param py
     *
     * @return
     */
    std::array<double, 3> ProjectionMatrix::directionToPosition(const double pi,
                                                                const double pj) const
    {
        Matrix toinvert = this->colSubMatrix(3);
        LUDoolittleForm lf = LUDoolittleForm::LUDecomposeDoolittle(toinvert, 1e-10);
        SquareMatrix inv = lf.inverseMatrix();
        Matrix P(3, 1, { pi, pj, 1.0 });
        Matrix v = inv * P;
        std::array<double, 3> x;
        x[0] = v(0, 0);
        x[1] = v(1, 0);
        x[2] = v(2, 0);
        // Want projection on VN to be 1
        std::array<double, 3> VN = directionVectorVN();
        double product = vectorDotProduct<3>(x, VN);
        return multiplyVectorByConstant<3>(x, 1.0 / product);
    }

    /**
     *Actual projection
     */

    void ProjectionMatrix::project0(
        const double X0, const double X1, const double X2, double* pi, double* pj) const
    {
        // We are using vector (X,Y,Z,0)+S and use fact that rows are orthogonal to S
        // So its sufficient to use just first three items
        double divideByMe = X0 * (*this)(2, 0) + X1 * (*this)(2, 1) + X2 * (*this)(2, 2);
        if(std::abs(divideByMe) < zeroPrecisionTolerance)
        {
            *pi = NAN;
            *pj = NAN;
        } else
        {
            *pi = (X0 * (*this)(0, 0) + X1 * (*this)(0, 1) + X2 * (*this)(0, 2)) / divideByMe;
            *pj = (X0 * (*this)(1, 0) + X1 * (*this)(1, 1) + X2 * (*this)(1, 2)) / divideByMe;
        }
    }

    void ProjectionMatrix::project0(
        const float X0, const float X1, const float X2, float* pi, float* pj) const
    {
        // We are using vector (X,Y,Z,0)+S and use fact that rows are orthogonal to S
        // So its sufficient to use just first three items
        double pxd, pyd;
        project0(double(X0), double(X1), double(X2), &pxd, &pyd);
        *pi = float(pxd);
        *pj = float(pyd);
    }

    void ProjectionMatrix::project(
        const double x0, const double x1, const double x2, double* pi, double* pj) const
    {
        double divideByMe
            = x0 * (*this)(2, 0) + x1 * (*this)(2, 1) + x2 * (*this)(2, 2) + (*this)(2, 3);
        if(std::abs(divideByMe) < zeroPrecisionTolerance)
        {
            *pi = NAN;
            *pj = NAN;
        } else
        {
            *pi = (x0 * (*this)(0, 0) + x1 * (*this)(0, 1) + x2 * (*this)(0, 2) + (*this)(0, 3))
                / divideByMe;
            *pj = (x0 * (*this)(1, 0) + x1 * (*this)(1, 1) + x2 * (*this)(1, 2) + (*this)(1, 3))
                / divideByMe;
        }
    }

    /**
     *Actual projection
     */
    void ProjectionMatrix::project(
        const float x0, const float x1, const float x2, float* pi, float* pj) const
    {
        double pid, pjd;
        project(double(x0), double(x1), double(x2), &pid, &pjd);
        *pi = float(pid);
        *pj = float(pjd);
    }

    ProjectionMatrix ProjectionMatrix::shiftDetectorOrigin(double x, double y) const
    {
        Matrix T = Matrix::unitDiagonal(3, 3);
        T(0, 2) = -x;
        T(1, 2) = -y;
        return T * (*this);
    }

    std::string ProjectionMatrix::toString(std::string name) const
    {
        std::ostringstream os;
        matrix::RQFactorization rq;
        std::shared_ptr<CTL::matrix::Matrix> F
            = std::make_shared<CTL::matrix::Matrix>(this->colSubMatrix(3));
        Matrix p4(3, 1, { (*this)(0, 3), (*this)(1, 3), (*this)(2, 3) });
        rq.factorize(F);
        auto C = rq.getRMatrix();
        auto Q = rq.getQMatrix();
        // I use backward substitution from the LUDoolittleForm, where C matrix is upper diagonal
        // representing LU matrix
        std::shared_ptr<std::vector<uint32_t>> P = std::make_shared<std::vector<uint32_t>>();
        for(int i = 0; i != 3; i++)
        {
            P->push_back(i);
            // On i-th position is the row index in the original matrix that is on the i-th
            // position in the PA matrix.
        }
        LUDoolittleForm lu(std::make_shared<SquareMatrix>(*C), P, false);
        Matrix u = lu.backwardSubstitute(p4);
        auto Qt = (-1.0) * Q->T();
        auto Qu = Qt * u;
        double factor = 1 / (*C)(2, 2);
        Matrix Cprime = factor * (*C);
        std::array<double, 3> S = sourcePosition();
        for(unsigned int i = 0; i != 3; ++i)
        {
            if(i == 1)
            {
                os << io::xprintf("%s = |", name.c_str());
            } else
            {
                os << "    |";
            }
            for(unsigned int j = 0; j != 4; ++j)
            {
                if(j != 0)
                    os << " ";
                os << std::setw(9) << std::fixed << std::setfill(' ') << std::setprecision(3)
                   << static_cast<double>((*this)(i, j));
            }
            os << "|";
            if(i == 1)
            {
                os << " = C[Q|u] = ";
            } else
            {
                os << "            ";
            }
            if(i == 1)
            {
                os << io::xprintf("%9.1f", factor);
            } else
            {
                os << "         ";
            }
            os << "|";
            for(unsigned int j = 0; j != 3; ++j)
            {
                if(j != 0)
                    os << " ";
                os << std::setw(9) << std::fixed << std::setfill(' ') << std::setprecision(3)
                   << static_cast<double>((Cprime)(i, j));
            }
            if(i == 1)
            {
                os << "|.|";

            } else
            {
                os << "| |";
            }
            for(unsigned int j = 0; j != 3; ++j)
            {
                if(j != 0)
                    os << " ";
                os << std::setw(9) << std::fixed << std::setfill(' ') << std::setprecision(3)
                   << static_cast<double>((*Q)(i, j));
            }
            os << "|";
            os << std::setw(9) << std::fixed << std::setfill(' ') << std::setprecision(3)
               << static_cast<double>(u(i, 0));
            os << "|";
            os << "\n";
        }
        os << io::xprintf("S = [%5.2f, %5.2f, %5.2f]", S[0], S[1], S[2]);
        os << io::xprintf(", -Q^T u = [%5.2f, %5.2f, %5.2f].", Qu(0, 0), Qu(1, 0), Qu(2, 0));
        os << std::endl;
        return os.str();
    }

    void ProjectionMatrix::computeSourcePosition() // Called during construction
    {
        std::array<double, 4> s;
        s[0] = LUDoolittleForm::determinant(this->colSubMatrix(0));
        s[1] = -LUDoolittleForm::determinant(this->colSubMatrix(1));
        s[2] = LUDoolittleForm::determinant(this->colSubMatrix(2));
        s[3] = -LUDoolittleForm::determinant(this->colSubMatrix(3));
        std::array<double, 4> r3 = { A[8], A[9], A[10], A[11] };
        s = matrix::orthogonalPartOfVectorWithRespectToSecondVector(s,
                                                                    r3); // Probably not effective
        source[0] = s[0] / s[3];
        source[1] = s[1] / s[3];
        source[2] = s[2] / s[3];
    }

    std::array<double, 2> ProjectionMatrix::principalRayProjection() const
    {
        std::array<double, 3> VN = directionVectorVN();
        std::array<double, 2> PRP;
        project0(VN[0], VN[1], VN[2], &PRP[0], &PRP[1]);
        return PRP;
    }

    std::array<double, 2> ProjectionMatrix::focalLength() const
    {
        std::array<double, 2> FL;
        FL[0] = matrix::vectorNorm(directionVectorVX());
        FL[1] = matrix::vectorNorm(directionVectorVY());
        return FL;
    }

    std::array<double, 2> ProjectionMatrix::pixelSizes(const double sourceToDetector) const
    {
        std::array<double, 2> FL = focalLength();
        std::array<double, 2> X;
        X[0] = sourceToDetector / FL[0];
        X[1] = sourceToDetector / FL[1];
        return X;
    }

    double ProjectionMatrix::sourceToDetectorFromPX(const double PX) const
    {
        return PX * matrix::vectorNorm(directionVectorVX());
    }

    double ProjectionMatrix::sourceToDetectorFromPY(const double PY) const
    {
        return PY * matrix::vectorNorm(directionVectorVY());
    }

    void ProjectionMatrix::sourcePosition(double* vector3) const
    {
        std::copy(std::begin(source), std::begin(source) + 3, vector3);
    }
    void ProjectionMatrix::normalToDetector(double* vector3) const

    {
        std::array<double, 3> V = normalToDetector();
        std::copy(std::begin(V), std::begin(V) + 3, vector3);
    }
    void ProjectionMatrix::principalRayProjection(double* vector2) const
    {
        std::array<double, 2> V = principalRayProjection();
        std::copy(std::begin(V), std::begin(V) + 2, vector2);
    }
    void ProjectionMatrix::focalLength(double* vector2) const
    {
        std::array<double, 2> V = focalLength();
        std::copy(std::begin(V), std::begin(V) + 2, vector2);
    }
    void ProjectionMatrix::pixelSizes(const double sourceToDetector, double* vector2) const
    {
        std::array<double, 2> V = pixelSizes(sourceToDetector);
        std::copy(std::begin(V), std::begin(V) + 2, vector2);
    }
    void
    ProjectionMatrix::directionToPosition(const double pi, const double pj, double* vector3) const
    {
        std::array<double, 3> V = directionToPosition(pi, pj);
        std::copy(std::begin(V), std::begin(V) + 3, vector3);
    }
    void ProjectionMatrix::projectionMatrixAsVector12(double* vector12) const
    {
        for(uint32_t i = 0; i != 3; i++)
        {
            for(uint32_t j = 0; j != 4; j++)
            {
                vector12[i * 4 + j] = (*this)(i, j);
            }
        }
    }

    void ProjectionMatrix::projectionMatrixAsVector9(double* vector9) const
    {
        for(uint32_t i = 0; i != 3; i++)
        {
            for(uint32_t j = 0; j != 3; j++)
            {
                vector9[i * 3 + j] = (*this)(i, j);
            }
        }
    }

    void ProjectionMatrix::inverseProjectionMatrixAsVector16(double* vector16) const
    {
        Matrix extendedMatrix(4, 4,
                              { (*this)(0, 0), (*this)(0, 1), (*this)(0, 2), (*this)(0, 3),
                                (*this)(1, 0), (*this)(1, 1), (*this)(1, 2), (*this)(1, 3),
                                (*this)(2, 0), (*this)(2, 1), (*this)(2, 2), (*this)(2, 3),
                                source[0], source[1], source[2], 1.0 });
        std::shared_ptr<CTL::matrix::SquareMatrix> ex
            = std::make_shared<CTL::matrix::SquareMatrix>(extendedMatrix);
        LUDoolittleForm lf = LUDoolittleForm::LUDecomposeDoolittle(*ex, 1e-10);
        SquareMatrix inv = lf.inverseMatrix();
        for(uint32_t i = 0; i != 4; i++)
        {
            for(uint32_t j = 0; j != 4; j++)
            {
                vector16[i * 4 + j] = inv(i, j);
            }
        }
    }

    void ProjectionMatrix::inverseProjectionMatrixAsVector9(double* vector9) const
    {
        SquareMatrix toinvert = this->colSubMatrix(3);
        LUDoolittleForm lf = LUDoolittleForm::LUDecomposeDoolittle(toinvert, 1e-10);
        SquareMatrix inv = lf.inverseMatrix();
        for(uint32_t i = 0; i != 3; i++)
        {
            for(uint32_t j = 0; j != 3; j++)
            {
                vector9[i * 3 + j] = inv(i, j);
            }
        }
    }
    std::array<double, 12> ProjectionMatrix::projectionMatrixAsVector12() const
    {
        std::array<double, 12> v;
        projectionMatrixAsVector12(std::begin(v));
        return v;
    }
    std::array<double, 9> ProjectionMatrix::projectionMatrixAsVector9() const
    {

        std::array<double, 9> v;
        projectionMatrixAsVector9(std::begin(v));
        return v;
    }
    std::array<double, 16> ProjectionMatrix::inverseProjectionMatrixAsVector16() const
    {

        std::array<double, 16> v;
        inverseProjectionMatrixAsVector16(std::begin(v));
        return v;
    }
    std::array<double, 9> ProjectionMatrix::inverseProjectionMatrixAsVector9() const
    {
        std::array<double, 9> v;
        inverseProjectionMatrixAsVector9(std::begin(v));
        return v;
    }

} // namespace matrix
} // namespace CTL
