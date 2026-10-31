#pragma once

#include "optpricing/types.hpp"

#include <cstddef>
#include <cstdint>

namespace optpricing {

struct MonteCarloConfig {
    std::size_t paths = 100000;   // requested number of simulated paths
    std::size_t steps = 252;      // barrier monitoring steps (European options use the exact terminal draw)
    unsigned threads = 1;         // worker threads; 0 = hardware concurrency
    std::uint64_t seed = 20260101;
    bool antithetic = true;       // pair each draw z with -z
    bool control_variate = true;  // discounted terminal price (European) / vanilla payoff (barrier)
};

struct McEstimate {
    double price;       // point estimate
    double std_error;   // standard error of the estimate
    std::size_t paths;  // number of simulated paths actually used
};

// European option by Monte Carlo under geometric Brownian motion.
// Reproducible for a fixed (seed, threads, paths, antithetic) configuration.
McEstimate mc_european(const Market& m, const VanillaOption& o, double vol, const MonteCarloConfig& cfg);

// Barrier option by Monte Carlo with discrete monitoring at `cfg.steps` dates. The barrier is
// shifted by the Broadie-Glasserman-Kou continuity correction, so the result approximates
// continuous monitoring. The residual bias is O(1/steps) and is covered by the tests.
McEstimate mc_barrier(const Market& m, const BarrierOption& o, double vol, const MonteCarloConfig& cfg);

}  // namespace optpricing
