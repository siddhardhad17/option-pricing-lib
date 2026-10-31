#include "optpricing/barrier.hpp"

#include "optpricing/black_scholes.hpp"
#include "optpricing/normal.hpp"

#include <algorithm>
#include <cmath>

namespace optpricing {
namespace {

// Building blocks of the Reiner-Rubinstein formulas.
struct Terms {
    double A;
    double B;
    double C;
    double D;
};

// phi = +1 for calls, -1 for puts; eta = +1 for down barriers, -1 for up barriers.
Terms barrier_terms(const Market& m, double strike, double barrier, double T, double vol,
                    double phi, double eta) {
    const double b = m.rate - m.dividend;
    const double s = vol * std::sqrt(T);
    const double mu = (b - 0.5 * vol * vol) / (vol * vol);

    const double x1 = std::log(m.spot / strike) / s + (1.0 + mu) * s;
    const double x2 = std::log(m.spot / barrier) / s + (1.0 + mu) * s;
    const double y1 = std::log(barrier * barrier / (m.spot * strike)) / s + (1.0 + mu) * s;
    const double y2 = std::log(barrier / m.spot) / s + (1.0 + mu) * s;

    const double spot_term = m.spot * std::exp((b - m.rate) * T);
    const double strike_term = strike * std::exp(-m.rate * T);
    const double ratio = barrier / m.spot;
    const double ratio_hi = std::pow(ratio, 2.0 * (mu + 1.0));
    const double ratio_lo = std::pow(ratio, 2.0 * mu);

    Terms t{};
    t.A = phi * spot_term * norm_cdf(phi * x1) - phi * strike_term * norm_cdf(phi * (x1 - s));
    t.B = phi * spot_term * norm_cdf(phi * x2) - phi * strike_term * norm_cdf(phi * (x2 - s));
    t.C = phi * spot_term * ratio_hi * norm_cdf(eta * y1) -
          phi * strike_term * ratio_lo * norm_cdf(eta * (y1 - s));
    t.D = phi * spot_term * ratio_hi * norm_cdf(eta * y2) -
          phi * strike_term * ratio_lo * norm_cdf(eta * (y2 - s));
    return t;
}

}  // namespace

double barrier_closed_form(const Market& m, const BarrierOption& o, double vol) {
    validate_market(m);
    validate_barrier(o);
    validate_vol(vol);

    const bool down = is_down_barrier(o.barrier);
    const bool knock_out = is_knock_out(o.barrier);
    const bool is_call = o.type == OptionType::Call;
    const double S = m.spot;
    const double H = o.barrier_level;
    const double K = o.strike;
    const double T = o.maturity;

    // Barrier already breached at inception.
    const bool breached = down ? (S <= H) : (S >= H);
    if (breached) {
        return knock_out ? 0.0 : bs_price(m, VanillaOption{o.type, K, T}, vol);
    }

    const double eta = down ? 1.0 : -1.0;
    const Terms c = barrier_terms(m, K, H, T, vol, +1.0, eta);
    const Terms p = barrier_terms(m, K, H, T, vol, -1.0, eta);
    const bool strike_above_barrier = K > H;

    double value = 0.0;
    if (down) {
        if (is_call) {
            if (knock_out) {
                value = strike_above_barrier ? c.A - c.C : c.B - c.D;
            } else {
                value = strike_above_barrier ? c.C : c.A - c.B + c.D;
            }
        } else {
            if (knock_out) {
                value = strike_above_barrier ? p.A - p.B + p.C - p.D : 0.0;
            } else {
                value = strike_above_barrier ? p.B - p.C + p.D : p.A;
            }
        }
    } else {
        if (is_call) {
            if (knock_out) {
                value = strike_above_barrier ? 0.0 : c.A - c.B + c.C - c.D;
            } else {
                value = strike_above_barrier ? c.A : c.B - c.C + c.D;
            }
        } else {
            if (knock_out) {
                value = strike_above_barrier ? p.B - p.D : p.A - p.C;
            } else {
                value = strike_above_barrier ? p.A - p.B + p.D : p.C;
            }
        }
    }
    // Removes rounding noise of order 1e-15 only; genuine negatives indicate a bug.
    return std::max(value, 0.0);
}

}  // namespace optpricing
