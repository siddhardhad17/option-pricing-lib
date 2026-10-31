#include "optpricing/normal.hpp"

#include <cmath>

namespace optpricing {
namespace {

constexpr double kInvSqrt2Pi = 0.398942280401432677939946059934;
constexpr double kInvSqrt2 = 0.707106781186547524400844362105;

}  // namespace

double norm_pdf(double x) noexcept {
    return kInvSqrt2Pi * std::exp(-0.5 * x * x);
}

double norm_cdf(double x) noexcept {
    return 0.5 * std::erfc(-x * kInvSqrt2);
}

}  // namespace optpricing
