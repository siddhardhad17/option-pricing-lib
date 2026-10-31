#pragma once

#include "optpricing/svi.hpp"

#include <limits>
#include <string>
#include <vector>

namespace optpricing {

// Log-moneyness grid on which the static checks are evaluated.
struct CheckGrid {
    double k_min = -2.0;
    double k_max = 2.0;
    int points = 801;
};

struct ArbitrageReport {
    bool passed = true;
    bool parameters_valid = true;     // every slice satisfies svi_params_valid
    bool positive_variance = true;    // w(k) > 0 on the grid
    bool lee_wing_bounds = true;      // b(1 + |rho|) <= 2 (Roger Lee moment bound)
    bool butterfly = true;            // Durrleman density g(k) >= 0 (no negative butterfly)
    bool calendar = true;             // w(T2, k) >= w(T1, k) for consecutive maturities
    double min_butterfly_density = std::numeric_limits<double>::infinity();
    double worst_calendar_gap = std::numeric_limits<double>::infinity();
    std::vector<std::string> violations;  // human-readable findings
};

// Durrleman's butterfly density g(k); the call-price density is non-negative iff g(k) >= 0.
double butterfly_density(const SviParams& p, double k) noexcept;

// Checks a surface given as slices sorted by strictly increasing maturity.
ArbitrageReport check_static_arbitrage(const std::vector<SviSlice>& slices, const CheckGrid& grid = {});

}  // namespace optpricing
