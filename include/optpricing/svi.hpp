#pragma once

#include "optpricing/types.hpp"

#include <vector>

namespace optpricing {

// Raw SVI parameterisation (Gatheral, 2004) of the total implied variance w = sigma_imp^2 * T
// as a function of log-moneyness k = ln(K / F):
//
//     w(k) = a + b * ( rho * (k - m) + sqrt((k - m)^2 + sigma^2) )
struct SviParams {
    double a;      // vertical level
    double b;      // angle between the wings (b >= 0)
    double rho;    // skew, |rho| < 1
    double m;      // horizontal translation
    double sigma;  // smoothness of the ATM vertex (sigma > 0)
};

struct SviSlice {
    double maturity;  // years, strictly increasing across a surface
    SviParams params;
};

// One market quote: a strike and the implied volatility observed for it.
struct SmileQuote {
    double strike;
    double implied_vol;
};

struct SviFitOptions {
    int max_iterations = 20000;
    int restarts = 3;            // extra Nelder-Mead passes from the incumbent optimum
    double tolerance = 1e-16;    // objective tolerance
};

struct SviFit {
    SviParams params;
    double rmse_total_variance;  // root-mean-square error in total variance
    double rmse_implied_vol;     // root-mean-square error in implied volatility
    int iterations;
    bool converged;              // optimiser converged and parameters satisfy the constraints
};

double svi_total_variance(const SviParams& p, double k) noexcept;
double svi_total_variance_slope(const SviParams& p, double k) noexcept;      // dw/dk
double svi_total_variance_curvature(const SviParams& p, double k) noexcept;  // d2w/dk2

// True if b >= 0, sigma > 0, |rho| < 1 and the minimum total variance is non-negative.
bool svi_params_valid(const SviParams& p) noexcept;

// Converts a market option price to a quote by inverting Black-Scholes.
SmileQuote quote_from_price(const Market& m, const VanillaOption& o, double price);

// Least-squares fit of one SVI slice to quotes. `forward` is F = S * exp((r - q) T).
// Needs at least five quotes with positive strikes and volatilities.
SviFit fit_svi_slice(const std::vector<SmileQuote>& quotes, double forward, double maturity,
                     const SviFitOptions& opt = {});

}  // namespace optpricing
