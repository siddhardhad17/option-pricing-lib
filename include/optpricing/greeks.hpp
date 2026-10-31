#pragma once

#include "optpricing/types.hpp"

#include <functional>

namespace optpricing {

// Any pricer as a function of (market, volatility, maturity).
using PriceFn = std::function<double(const Market& market, double vol, double maturity)>;

// Central finite-difference Greeks for an arbitrary pricer. Step sizes are scaled to each input.
Greeks finite_difference_greeks(const PriceFn& price, const Market& market, double vol, double maturity);

Greeks fd_greeks_vanilla(const Market& market, const VanillaOption& option, double vol);
Greeks fd_greeks_barrier(const Market& market, const BarrierOption& option, double vol);

}  // namespace optpricing
