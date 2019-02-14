#pragma once
// Internal
#include "MATRIX/ProjectionMatrix.hpp"
#include "MATRIX/ProjectionMatrixComparatorI.hpp"

namespace CTL {
namespace matrix {
    /**
    Interface to compare two projection matrices.
    */
    class SourcePositionComparator : public ProjectionMatrixComparatorI
    {
    public:
        double distance(ProjectionMatrix a, ProjectionMatrix b) override;
        /**Get the distance between two compared matrices*/
    };
} // namespace matrix
} // namespace CTL
