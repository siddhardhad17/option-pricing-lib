// Builds a synthetic three-expiry market, prices the quotes, inverts them to implied vols, fits
// an SVI slice per expiry, and runs the static no-arbitrage checks on the fitted surface.

#include "optpricing/black_scholes.hpp"
#include "optpricing/no_arbitrage.hpp"
#include "optpricing/svi.hpp"

#include <cmath>
#include <exception>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace optpricing;

int main() {
    try {
        const Market market{100.0, 0.03, 0.01};

        // "True" surface used to generate the synthetic market. Total variance grows with
        // maturity at every strike (the level a rises, the shape is constant), so it is
        // calendar-consistent by construction.
        const std::vector<SviSlice> truth = {
            {0.25, {0.012, 0.30, -0.40, 0.00, 0.18}},
            {0.50, {0.026, 0.30, -0.40, 0.00, 0.18}},
            {1.00, {0.054, 0.30, -0.40, 0.00, 0.18}},
        };

        std::cout << std::fixed << std::setprecision(6);
        std::vector<SviSlice> fitted;
        for (const SviSlice& slice : truth) {
            const double T = slice.maturity;
            const double forward = market.spot * std::exp((market.rate - market.dividend) * T);

            // Quotes: option prices on a strike grid, converted to implied vols.
            std::vector<SmileQuote> quotes;
            for (int i = -6; i <= 6; ++i) {
                const double k = 0.05 * i;
                const double K = forward * std::exp(k);
                const VanillaOption call{OptionType::Call, K, T};
                const double vol = std::sqrt(svi_total_variance(slice.params, k) / T);
                quotes.push_back(quote_from_price(market, call, bs_price(market, call, vol)));
            }

            const SviFit fit = fit_svi_slice(quotes, forward, T);
            std::cout << "T=" << T << "  converged=" << (fit.converged ? "yes" : "no")
                      << "  iterations=" << fit.iterations << "\n"
                      << "  params  a=" << fit.params.a << " b=" << fit.params.b << " rho=" << fit.params.rho
                      << " m=" << fit.params.m << " sigma=" << fit.params.sigma << "\n"
                      << "  rmse    total variance=" << std::scientific << std::setprecision(2)
                      << fit.rmse_total_variance << "  implied vol=" << fit.rmse_implied_vol << std::fixed
                      << std::setprecision(6) << "\n";
            fitted.push_back({T, fit.params});
        }

        const ArbitrageReport report = check_static_arbitrage(fitted);
        std::cout << "\nStatic no-arbitrage checks: " << (report.passed ? "PASSED" : "FAILED") << "\n"
                  << "  min butterfly density g(k) = " << report.min_butterfly_density << "\n"
                  << "  worst calendar gap         = " << report.worst_calendar_gap << "\n";
        for (const auto& v : report.violations) {
            std::cout << "  - " << v << "\n";
        }
        return report.passed ? 0 : 2;
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << "\n";
        return 1;
    }
}
