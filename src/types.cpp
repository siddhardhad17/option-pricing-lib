#include "optpricing/types.hpp"

#include <cmath>
#include <stdexcept>

namespace optpricing {
namespace {

void require(bool ok, const char* message) {
    if (!ok) {
        throw std::invalid_argument(message);
    }
}

bool positive_finite(double x) noexcept {
    return std::isfinite(x) && x > 0.0;
}

}  // namespace

void validate_market(const Market& m) {
    require(positive_finite(m.spot), "spot must be finite and > 0");
    require(std::isfinite(m.rate), "rate must be finite");
    require(std::isfinite(m.dividend), "dividend must be finite");
}

void validate_vol(double vol) {
    require(positive_finite(vol), "volatility must be finite and > 0");
}

void validate_vanilla(const VanillaOption& o) {
    require(positive_finite(o.strike), "strike must be finite and > 0");
    require(positive_finite(o.maturity), "maturity must be finite and > 0");
}

void validate_barrier(const BarrierOption& o) {
    require(positive_finite(o.strike), "strike must be finite and > 0");
    require(positive_finite(o.barrier_level), "barrier level must be finite and > 0");
    require(positive_finite(o.maturity), "maturity must be finite and > 0");
}

}  // namespace optpricing
