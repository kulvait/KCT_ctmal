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

    std::array<double, 3> ProjectionMatrix::tangentToDetectorXDirection() const
    {
        std::array<double, 3> tangent;
        std::array<double, 4> thirdRow = { A[8], A[9], A[10], A[11] };
        std::array<double, 4> firstRow = { A[0], A[1], A[2], A[3] };
        std::array<double, 4> tangentF = reorthogonalize(firstRow, thirdRow);
        if(std::abs(tangentF[3]) < zeroPrecisionTolerance)
        { // Vector that projects to infinity is from infinity itself
          // If the fourth coordinate is zero, then I can add the 3 vector without fourth position
          // to source location and obtain finite tangent
            tangent[0] = tangentF[0];
            tangent[1] = tangentF[1];
            tangent[2] = tangentF[2];
        } else
        {
            tangent[0] = tangentF[0] / tangentF[3];
            tangent[1] = tangentF[1] / tangentF[3];
            tangent[2] = tangentF[2] / tangentF[3];
            tangent[0] = tangent[0] - source[0];
            tangent[1] = tangent[1] - source[1];
            tangent[2] = tangent[2] - source[2];
        }
        return normalizeVector(tangent);
    }

    std::array<double, 3> ProjectionMatrix::tangentToDetectorYDirection() const
    {
        std::array<double, 3> tangent;
        std::array<double, 4> thirdRow = { A[8], A[9], A[10], A[11] };
        std::array<double, 4> secondRow = { A[4], A[5], A[6], A[7] };
        std::array<double, 4> tangentF = reorthogonalize(secondRow, thirdRow);
        if(std::abs(tangentF[3]) < zeroPrecisionTolerance)
        { // Vector that projects to infinity is from infinity itself
          // If the fourth coordinate is zero, then I can add the 3 vector without fourth position
          // to source location and obtain finite tangent
            tangent[0] = tangentF[0];
            tangent[1] = tangentF[1];
            tangent[2] = tangentF[2];
        } else
        {
            tangent[0] = tangentF[0] / tangentF[3];
            tangent[1] = tangentF[1] / tangentF[3];
            tangent[2] = tangentF[2] / tangentF[3];
            // I have found a vector in the 0 based coordinates that projects to infinity
            tangent[0] = tangent[0] - source[0];
            tangent[1] = tangent[1] - source[1];
            tangent[2] = tangent[2] - source[2];
        }
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
        /*Previous implementation
        std::array<double, 3> normal;
        normal[0] = A[8]; // (*this)(2, 0);
        normal[1] = A[9]; //(*this)(2, 1);
        normal[2] = A[10]; //(*this)(2, 2);
        */
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
    std::array<double, 3> ProjectionMatrix::projectedToPosition(double px, double py) const
    {
        Matrix rot = this->colSubMatrix(3);
        Matrix shift = Matrix(3, 3, { 1.0, 0.0, -px, 0.0, 1.0, -py, 0.0, 0.0, 1.0 });
        Matrix B = shift * rot; // First part of projection matrix that projects (x,y,z) to (0,0,?)
        std::array<double, 3> normal;
        // normal[0] = LUDoolittleForm::determinant(B.minorSubMatrix(2, 0));
        // normal[1] = -LUDoolittleForm::determinant(B.minorSubMatrix(2, 1));
        // normal[2] = LUDoolittleForm::determinant(B.minorSubMatrix(2, 2));
        normal[0] = B(0, 1) * B(1, 2) - B(0, 2) * B(1, 1);
        normal[1] = B(0, 2) * B(1, 0) - B(0, 0) * B(1, 2);
        normal[2] = B(0, 0) * B(1, 1) - B(0, 1) * B(1, 0);
        // Renormalization
        double normvec
            = std::sqrt(normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]);
        normal[0] /= normvec;
        normal[1] /= normvec;
        normal[2] /= normvec;
        return normal;
    }

    /**
     *Actual projection
     */
    void ProjectionMatrix::project(double x, double y, double z, double* px, double* py)
    {
        double divideByMe
            = x * (*this)(2, 0) + y * (*this)(2, 1) + z * (*this)(2, 2) + (*this)(2, 3);
        if(std::abs(divideByMe) < zeroPrecisionTolerance)
        {
            *px = NAN;
            *py = NAN;
        } else
        {
            *px = (x * (*this)(0, 0) + y * (*this)(0, 1) + z * (*this)(0, 2) + (*this)(0, 3))
                / divideByMe;
            *py = (x * (*this)(1, 0) + y * (*this)(1, 1) + z * (*this)(1, 2) + (*this)(1, 3))
                / divideByMe;
        }
    }

    /**
     *Actual projection
     */
    void ProjectionMatrix::project(float x, float y, float z, float* px, float* py)
    {
        double divideByMe
            = x * (*this)(2, 0) + y * (*this)(2, 1) + z * (*this)(2, 2) + (*this)(2, 3);
        if(std::abs(divideByMe) < zeroPrecisionTolerance)
        {
            *px = NAN;
            *py = NAN;
        } else
        {
            double pxd = (x * (*this)(0, 0) + y * (*this)(0, 1) + z * (*this)(0, 2) + (*this)(0, 3))
                / divideByMe;
            double pyd = (x * (*this)(1, 0) + y * (*this)(1, 1) + z * (*this)(1, 2) + (*this)(1, 3))
                / divideByMe;
            *px = float(pxd);
            *py = float(pyd);
        }
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

    double ProjectionMatrix::vectorNorm(std::array<double, 3> v) const
    {
        double normsq = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
        return std::sqrt(normsq);
    }

    std::array<double, 3> ProjectionMatrix::normalizeVector(std::array<double, 3> v) const
    {
        std::array<double, 3> normal;
        double normvec = vectorNorm(v);
        normal[0] = v[0] / normvec;
        normal[1] = v[1] / normvec;
        normal[2] = v[2] / normvec;
        return normal;
    }

    std::array<double, 4> ProjectionMatrix::reorthogonalize(std::array<double, 4> v,
                                                            std::array<double, 4> og) const
    {
        double nogsq = og[0] * og[0] + og[1] * og[1] + og[2] * og[2] + og[3] * og[3];
        double vog = v[0] * og[0] + v[1] * og[1] + v[2] * og[2] + v[3] * og[3];
        double factor = vog / nogsq;
        std::array<double, 4> orthogonalized;
        for(int i = 0; i != 4; i++)
        {
            orthogonalized[i] = v[i] - factor * og[i];
        }
        return orthogonalized;
    }

    void ProjectionMatrix::computeSourcePosition() // Called during construction
    {
        std::array<double, 4> s;
        s[0] = LUDoolittleForm::determinant(this->colSubMatrix(0));
        s[1] = -LUDoolittleForm::determinant(this->colSubMatrix(1));
        s[2] = LUDoolittleForm::determinant(this->colSubMatrix(2));
        s[3] = -LUDoolittleForm::determinant(this->colSubMatrix(3));
        std::array<double, 4> r3 = { A[8], A[9], A[10], A[11] };
        s = reorthogonalize(s, r3); // Probably not effective
        source[0] = s[0] / s[3];
        source[1] = s[1] / s[3];
        source[2] = s[2] / s[3];
    }

} // namespace matrix
} // namespace CTL
