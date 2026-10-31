#include "optpricing/greeks.hpp"

#include "optpricing/barrier.hpp"
#include "optpricing/black_scholes.hpp"

#include <stdexcept>

namespace optpricing {

Greeks finite_difference_greeks(const PriceFn& price, const Market& market, double vol, double maturity) {
    validate_market(market);
    validate_vol(vol);
    if (!(maturity > 1e-3)) {
        throw std::invalid_argument("maturity must be > 1e-3 for finite differences");
    }
    if (!(vol > 1e-3)) {
        throw std::invalid_argument("volatility must be > 1e-3 for finite differences");
    }

    const double S = market.spot;
    const double h_spot = 1e-4 * S;   // delta: truncation O(h^2), round-off ~ eps / h
    const double h_gamma = 1e-3 * S;  // gamma needs a larger step because it divides by h^2
    const double h_vol = 1e-4;
    const double h_time = 1e-4;
    const double h_rate = 1e-4;

    auto with_spot = [&](double s) {
        Market mm = market;
        mm.spot = s;
        return mm;
    };
    auto with_rate = [&](double r) {
        Market mm = market;
        mm.rate = r;
        return mm;
    };

    const double v0 = price(market, vol, maturity);
    Greeks g{};
    g.delta = (price(with_spot(S + h_spot), vol, maturity) - price(with_spot(S - h_spot), vol, maturity)) /
              (2.0 * h_spot);
    g.gamma = (price(with_spot(S + h_gamma), vol, maturity) - 2.0 * v0 +
               price(with_spot(S - h_gamma), vol, maturity)) /
              (h_gamma * h_gamma);
    g.vega = (price(market, vol + h_vol, maturity) - price(market, vol - h_vol, maturity)) / (2.0 * h_vol);
    // theta = -dV/dT: time decay is the negative of the derivative with respect to maturity.
    g.theta = -(price(market, vol, maturity + h_time) - price(market, vol, maturity - h_time)) / (2.0 * h_time);
    g.rho = (price(with_rate(market.rate + h_rate), vol, maturity) -
             price(with_rate(market.rate - h_rate), vol, maturity)) /
            (2.0 * h_rate);
    return g;
}

Greeks fd_greeks_vanilla(const Market& market, const VanillaOption& option, double vol) {
    return finite_difference_greeks(
        [option](const Market& m, double v, double t) {
            VanillaOption o = option;
            o.maturity = t;
            return bs_price(m, o, v);
        },
        market, vol, option.maturity);
}

Greeks fd_greeks_barrier(const Market& market, const BarrierOption& option, double vol) {
    return finite_difference_greeks(
        [option](const Market& m, double v, double t) {
            BarrierOption o = option;
            o.maturity = t;
            return barrier_closed_form(m, o, v);
        },
        market, vol, option.maturity);
}

}  // namespace optpricing
