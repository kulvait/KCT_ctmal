#pragma once

#include "MATRIX/ProjectionMatrix.hpp"
#include "ProjectionMatrixComparatorI.hpp" //utils.h
#include "matrix.h"
#include "pmatcomparator.h"

namespace CTL {
namespace matrix {

    class ProjectionErrorComparator : public ProjectionMatrixComparatorI
    {
    private:
        bool useMaxErrorInsteadOfMeanError;

    public:
        ProjectionErrorComparator(bool useMaxErrorInsteadOfMeanError)
        {
            this->useMaxErrorInsteadOfMeanError = useMaxErrorInsteadOfMeanError;
        }

        double distance(ProjectionMatrix a, ProjectionMatrix b) override
        {
            PMatComparator pmc;
CTL::Matrix<3,4> a_;
CTL::Matrix<3,4> b_;
for(int i = 0; i!=3 ; i++)
	for(int j=0; j!=4; j++)
	{
		a_(i,j) = a(i,j);
		b_(i,j) = b(i,j);
	}

            auto eval = pmc(a_, b_);
            if(useMaxErrorInsteadOfMeanError)
            {
                return eval.maxError;
            } else
            {
                return eval.meanError;
            }
        }
    };
} // namespace util
} // namespace CTL
