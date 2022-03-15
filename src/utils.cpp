#include "MATRIX/utils.hpp"

namespace KCT {
namespace matrix {

std::array<double, 3> vectorProduct(std::array<double, 3> v,
                                    std::array<double, 3> w) {
  std::array<double, 3> vp;
  vp[0] = v[1] * w[2] - v[2] * w[1];
  vp[1] = v[2] * w[0] - v[0] * w[2];
  vp[2] = v[0] * w[1] - v[1] * w[0];
  return vp;
}

} // namespace matrix
} // namespace KCT
