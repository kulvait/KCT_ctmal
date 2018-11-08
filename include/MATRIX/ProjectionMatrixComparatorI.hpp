#pragma once

#include "MATRIX/ProjectionMatrix.hpp"

namespace CTL {
namespace util {
    /**
    Interface to compare two projection matrices.
    */
    class ProjectionMatrixComparatorI
    {
    public:
        virtual double distance(ProjectionMatrix a, ProjectionMatrix b) = 0;
        /**Get the distance between two compared matrices*/
    };
} // namespace util
} // namespace CTL
