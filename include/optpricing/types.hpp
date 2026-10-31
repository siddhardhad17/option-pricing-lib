#pragma once

#include <cmath>
#include <cstddef>

namespace optpricing {

enum class OptionType { Call, Put };

enum class BarrierType { DownAndOut, DownAndIn, UpAndOut, UpAndIn };

// Market state shared by all pricers. Rates are continuously compounded and annualised.
struct Market {
    double spot;      // S0 > 0
    double rate;      // r, risk-free rate (any finite sign)
    double dividend;  // q, continuous dividend yield (any finite sign)
};

struct VanillaOption {
    OptionType type;
    double strike;    // K > 0
    double maturity;  // T > 0, in years
};

struct BarrierOption {
    OptionType type;
    BarrierType barrier;   // knock-in / knock-out, down / up
    double strike;         // K > 0
    double barrier_level;  // H > 0 (continuously monitored, zero rebate)
    double maturity;       // T > 0, in years
};

// Sensitivity conventions, shared by the analytic and finite-difference code:
//   delta = dV/dS,  gamma = d2V/dS2,  vega = dV/dsigma (per 1.00 of vol, not per 1%),
//   theta = dV/dt = -dV/dT (per year),  rho = dV/dr (per 1.00 of rate).
struct Greeks {
    double delta;
    double gamma;
    double vega;
    double theta;
    double rho;
};

// +1 for calls, -1 for puts.
inline double option_sign(OptionType t) noexcept {
    return t == OptionType::Call ? 1.0 : -1.0;
}

inline bool is_down_barrier(BarrierType b) noexcept {
    return b == BarrierType::DownAndOut || b == BarrierType::DownAndIn;
}

inline bool is_knock_out(BarrierType b) noexcept {
    return b == BarrierType::DownAndOut || b == BarrierType::UpAndOut;
}

// All validators throw std::invalid_argument on violated preconditions.
void validate_market(const Market& m);
void validate_vol(double vol);
void validate_vanilla(const VanillaOption& o);
void validate_barrier(const BarrierOption& o);

}  // namespace optpricing
