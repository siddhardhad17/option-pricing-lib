#include "optpricing/black_scholes.hpp"

#include "optpricing/normal.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace optpricing {
namespace {

struct D12 {
    double d1;
    double d2;
};

D12 compute_d(const Market& m, double strike, double maturity, double vol) {
    const double sd = vol * std::sqrt(maturity);
    const double d1 = (std::log(m.spot / strike) + (m.rate - m.dividend + 0.5 * vol * vol) * maturity) / sd;
    return {d1, d1 - sd};
}

}  // namespace

double bs_price(const Market& m, const VanillaOption& o, double vol) {
    validate_market(m);
    validate_vanilla(o);
    validate_vol(vol);

    const D12 d = compute_d(m, o.strike, o.maturity, vol);
    const double disc_r = std::exp(-m.rate * o.maturity);
    const double disc_q = std::exp(-m.dividend * o.maturity);

    if (o.type == OptionType::Call) {
        return m.spot * disc_q * norm_cdf(d.d1) - o.strike * disc_r * norm_cdf(d.d2);
    }
    return o.strike * disc_r * norm_cdf(-d.d2) - m.spot * disc_q * norm_cdf(-d.d1);
}

Greeks bs_greeks(const Market& m, const VanillaOption& o, double vol) {
    validate_market(m);
    validate_vanilla(o);
    validate_vol(vol);

    const double T = o.maturity;
    const double sqrt_t = std::sqrt(T);
    const double disc_r = std::exp(-m.rate * T);
    const double disc_q = std::exp(-m.dividend * T);
    const D12 d = compute_d(m, o.strike, T, vol);
    const double pdf1 = norm_pdf(d.d1);
    const double sgn = option_sign(o.type);

    Greeks g{};
    g.delta = sgn * disc_q * norm_cdf(sgn * d.d1);
    g.gamma = disc_q * pdf1 / (m.spot * vol * sqrt_t);
    g.vega = m.spot * disc_q * pdf1 * sqrt_t;

    // theta = -dV/dT. Written with the sign convention shared by calls (sgn = +1) and puts (sgn = -1).
    const double decay = -m.spot * disc_q * pdf1 * vol / (2.0 * sqrt_t);
    const double carry = m.rate * o.strike * disc_r * norm_cdf(sgn * d.d2) -
                         m.dividend * m.spot * disc_q * norm_cdf(sgn * d.d1);
    g.theta = decay - sgn * carry;

    g.rho = sgn * o.strike * T * disc_r * norm_cdf(sgn * d.d2);
    return g;
}

double bs_implied_vol(const Market& m, const VanillaOption& o, double price, double tol) {
    validate_market(m);
    validate_vanilla(o);
    if (!(tol > 0.0)) {
        throw std::invalid_argument("tolerance must be > 0");
    }

    const double forward_disc = m.spot * std::exp(-m.dividend * o.maturity);
    const double strike_disc = o.strike * std::exp(-m.rate * o.maturity);
    double lower = 0.0;
    double upper = 0.0;
    if (o.type == OptionType::Call) {
        lower = std::max(forward_disc - strike_disc, 0.0);
        upper = forward_disc;
    } else {
        lower = std::max(strike_disc - forward_disc, 0.0);
        upper = strike_disc;
    }
    if (!(price > lower && price < upper)) {
        throw std::invalid_argument("option price violates no-arbitrage bounds");
    }

    constexpr double kLoVol = 1e-8;
    constexpr double kHiVol = 5.0;
    if (bs_price(m, o, kHiVol) < price) {
        throw std::invalid_argument("option price requires volatility above 500%");
    }

    // The Black-Scholes price is strictly increasing in volatility, so bisection is safe.
    double lo = kLoVol;
    double hi = kHiVol;
    for (int i = 0; i < 300 && hi - lo > tol; ++i) {
        const double mid = 0.5 * (lo + hi);
        if (bs_price(m, o, mid) < price) {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    return 0.5 * (lo + hi);
}

}  // namespace optpricing
