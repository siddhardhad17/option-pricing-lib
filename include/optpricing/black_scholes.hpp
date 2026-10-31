#pragma once

#include "optpricing/types.hpp"

namespace optpricing {

// Black-Scholes-Merton price of a European option with continuous dividend yield.
double bs_price(const Market& m, const VanillaOption& o, double vol);

// Closed-form Greeks (see the sign conventions documented in types.hpp).
Greeks bs_greeks(const Market& m, const VanillaOption& o, double vol);

// Implied volatility by bisection on the monotone pricing function.
// Throws std::invalid_argument if `price` lies outside the no-arbitrage bounds.
double bs_implied_vol(const Market& m, const VanillaOption& o, double price, double tol = 1e-12);

}  // namespace optpricing
