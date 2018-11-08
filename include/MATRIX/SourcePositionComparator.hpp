#pragma once
// Internal
#include "MATRIX/ProjectionMatrix.hpp"
#include "utils/ProjectionMatrixComparatorI.hpp"

namespace CTL {
namespace util {
    /**
    Interface to compare two projection matrices.
    */
    class SourcePositionComparator : public ProjectionMatrixComparatorI
    {
    public:
        double distance(ProjectionMatrix a, ProjectionMatrix b) override;
        /**Get the distance between two compared matrices*/
    };
} // namespace util
} // namespace CTL
