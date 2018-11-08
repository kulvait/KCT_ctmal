#include "MATRIX/SourcePositionComparator.hpp"

namespace CTL {
namespace util {

    double SourcePositionComparator::distance(ProjectionMatrix a, ProjectionMatrix b)
    {
        std::array<double, 3> s1 = a.sourcePosition();
        std::array<double, 3> s2 = b.sourcePosition();
        double square = 0;
        for(int i = 0; i != 3; i++)
        {
            square += (s1[i] - s2[i]) * (s1[i] - s2[i]);
        }
        return square;
    }
} // namespace util
} // namespace CTL
