#pragma once

#include "MATRIX/ProjectionMatrix.hpp"
#include "ProjectionMatrixComparatorI.hpp" //utils.h
#include "pmatcomparator.h"

namespace CTL {
namespace util {

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
            auto eval = pmc(a, b);
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
