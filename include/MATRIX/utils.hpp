#pragma once

#include <array>
#include <cmath>

namespace CTL {
namespace matrix {

    /**
     * @return Dot product of two vectors
     */
    template <std::size_t N>
    double vectorDotProduct(std::array<double, N> v, std::array<double, N> w);
    /**
     * @return l2 norm of a vector v
     */
    template <std::size_t N>
    double vectorNorm(std::array<double, N> v);
    /**
     * @return v + w
     */
    template <std::size_t N>
    std::array<double, N> vectorSum(std::array<double, N> v, std::array<double, N> w);
    /**
     * @return v - w
     */
    template <std::size_t N>
    std::array<double, N> vectorDiff(std::array<double, N> v, std::array<double, N> w);
    /**
     * @return Unit vector with the direction of v
     */
    template <std::size_t N>
    std::array<double, N> normalizeVector(std::array<double, N> v);
    /**
     * @return Unit vector with the direction of v
     */
    template <std::size_t N>
    std::array<double, N> multiplyVectorByConstant(std::array<double, N> v, const double c);
    /**
     * @return Part of v orthogonal to og
     */
    template <std::size_t N>
    std::array<double, N> orthogonalPartOfVectorWithRespectToSecondVector(std::array<double, N> v,
                                                                          std::array<double, N> og);
    // Definitions in header as we have templates
    template <std::size_t N>
    double vectorDotProduct(std::array<double, N> v, std::array<double, N> w)
    {
        double product = 0.0;
        for(std::size_t i = 0; i != N; i++)
        {
            product += (v[i] * w[i]);
        }
        return product;
    }

    template <std::size_t N>
    double vectorNorm(std::array<double, N> v)
    {
        double normsq = 0.0;
        for(std::size_t i = 0; i != N; i++)
        {
            normsq += (v[i] * v[i]);
        }
        return std::sqrt(normsq);
    }

    template <std::size_t N>
    std::array<double, N> vectorSum(std::array<double, N> v, std::array<double, N> w)
    {
        std::array<double, N> sum;
        for(std::size_t i = 0; i != N; i++)
        {
            sum[i] = v[i] + w[i];
        }
        return sum;
    }

    template <std::size_t N>
    std::array<double, N> vectorDiff(std::array<double, N> v, std::array<double, N> w)
    {
        std::array<double, N> diff;
        for(std::size_t i = 0; i != N; i++)
        {
            diff[i] = v[i] - w[i];
        }
        return diff;
    }

    template <std::size_t N>
    std::array<double, N> normalizeVector(std::array<double, N> v)
    {
        std::array<double, N> unitVector;
        double n = vectorNorm<N>(v);
        for(std::size_t i = 0; i != N; i++)
        {
            unitVector[i] = v[i] / n;
        }
        return unitVector;
    }

    template <std::size_t N>
    std::array<double, N> multiplyVectorByConstant(std::array<double, N> v, const double c)
    {
        std::array<double, N> x;
        for(std::size_t i = 0; i != N; i++)
        {
            x[i] = c * v[i];
        }
        return x;
    }

    template <std::size_t N>
    std::array<double, N> orthogonalPartOfVectorWithRespectToSecondVector(std::array<double, N> v,
                                                                          std::array<double, N> og)
    {
        double nogsq = vectorDotProduct(og, og);
        double vog = vectorDotProduct(v, og);
        double factor = vog / nogsq;
        std::array<double, N> orthogonalized;
        for(std::size_t i = 0; i != N; i++)
        {
            orthogonalized[i] = v[i] - factor * og[i];
        }
        return orthogonalized;
    }
} // namespace matrix
} // namespace CTL
