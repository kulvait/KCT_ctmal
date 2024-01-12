#pragma once

#include "MATRIX/ProjectionMatrix.hpp"

namespace KCT {
namespace matrix {
    /**
    Interface to compare two projection matrices.
    */
    class ProjectionMatrixComparatorI
    {
    public:
        virtual double distance(ProjectionMatrix a, ProjectionMatrix b) = 0;
        /**Get the distance between two compared matrices*/
        virtual ~ProjectionMatrixComparatorI() = default;
    };
} // namespace matrix
} // namespace KCT
