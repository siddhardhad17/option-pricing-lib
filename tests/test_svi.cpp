#include "test_harness.hpp"

#include "optpricing/black_scholes.hpp"
#include "optpricing/nelder_mead.hpp"
#include "optpricing/no_arbitrage.hpp"
#include "optpricing/svi.hpp"

#include <cmath>
#include <vector>

using namespace optpricing;

TEST_CASE(svi_formula_identities) {
    const SviParams p{0.04, 0.4, -0.4, 0.1, 0.3};
    // At k = m the total variance is a + b * sigma.
    CHECK_NEAR(svi_total_variance(p, p.m), p.a + p.b * p.sigma, 0.0, 1e-14);
    // Analytic slope and curvature agree with central differences.
    const double h = 1e-5;
    for (double k : {-0.5, 0.0, 0.3}) {
        const double d1 = (svi_total_variance(p, k + h) - svi_total_variance(p, k - h)) / (2 * h);
        const double d2 = (svi_total_variance(p, k + h) - 2 * svi_total_variance(p, k) + svi_total_variance(p, k - h)) /
                          (h * h);
        CHECK_NEAR(svi_total_variance_slope(p, k), d1, 1e-6, 1e-7);
        CHECK_NEAR(svi_total_variance_curvature(p, k), d2, 1e-4, 1e-5);
    }
    CHECK(svi_params_valid(p));
    CHECK(!svi_params_valid(SviParams{0.04, -0.1, 0.0, 0.0, 0.3}));  // negative b
    CHECK(!svi_params_valid(SviParams{0.04, 0.4, 1.0, 0.0, 0.3}));   // |rho| = 1
}

TEST_CASE(svi_fit_recovers_synthetic_slice) {
    // Build quotes from known SVI parameters, round-trip through prices, then fit.
    const Market m{100.0, 0.03, 0.01};
    const double T = 0.75;
    const double F = m.spot * std::exp((m.rate - m.dividend) * T);
    const SviParams truth{0.02, 0.35, -0.45, 0.05, 0.25};

    std::vector<SmileQuote> quotes;
    for (int i = -8; i <= 8; ++i) {
        const double k = 0.06 * i;
        const double K = F * std::exp(k);
        const double vol = std::sqrt(svi_total_variance(truth, k) / T);
        const double price = bs_price(m, {OptionType::Call, K, T}, vol);
        quotes.push_back(quote_from_price(m, {OptionType::Call, K, T}, price));
    }

    const SviFit fit = fit_svi_slice(quotes, F, T);
    CHECK(fit.converged);
    CHECK(fit.rmse_total_variance < 1e-6);
    CHECK(fit.rmse_implied_vol < 1e-5);
    CHECK(svi_params_valid(fit.params));
}

TEST_CASE(svi_fit_rejects_too_few_quotes) {
    const std::vector<SmileQuote> two = {{90.0, 0.2}, {110.0, 0.22}};
    CHECK_THROWS_INVALID(fit_svi_slice(two, 100.0, 1.0));
}

TEST_CASE(healthy_surface_passes_static_checks) {
    const std::vector<SviSlice> surface = {
        {0.25, {0.010, 0.30, -0.35, 0.00, 0.20}},
        {0.50, {0.025, 0.30, -0.35, 0.00, 0.20}},
        {1.00, {0.050, 0.30, -0.35, 0.00, 0.20}},
    };
    const ArbitrageReport rep = check_static_arbitrage(surface);
    CHECK(rep.passed);
    CHECK(rep.butterfly);
    CHECK(rep.calendar);
    CHECK(rep.lee_wing_bounds);
    CHECK(rep.violations.empty());
}

TEST_CASE(calendar_violation_is_detected) {
    // The longer slice has lower total variance at every strike.
    const std::vector<SviSlice> surface = {
        {0.5, {0.05, 0.30, -0.35, 0.00, 0.20}},
        {1.0, {0.02, 0.30, -0.35, 0.00, 0.20}},
    };
    const ArbitrageReport rep = check_static_arbitrage(surface);
    CHECK(!rep.passed);
    CHECK(!rep.calendar);
    CHECK(rep.worst_calendar_gap < 0.0);
    CHECK(!rep.violations.empty());
}

TEST_CASE(lee_bound_violation_is_detected) {
    // b(1 + rho) = 3.75 exceeds the asymptotic slope bound of 2.
    const std::vector<SviSlice> surface = {{1.0, {0.01, 2.5, 0.5, 0.0, 0.1}}};
    const ArbitrageReport rep = check_static_arbitrage(surface);
    CHECK(!rep.passed);
    CHECK(!rep.lee_wing_bounds);
}

TEST_CASE(butterfly_violation_is_detected) {
    // Satisfies the Lee bound (b(1 + |rho|) = 1.9 <= 2) but has a sharp vertex, so the density
    // g(k) turns negative near the money.
    const std::vector<SviSlice> surface = {{1.0, {0.01, 1.0, -0.9, 0.0, 0.02}}};
    const ArbitrageReport rep = check_static_arbitrage(surface);
    CHECK(rep.lee_wing_bounds);
    CHECK(!rep.butterfly);
    CHECK(rep.min_butterfly_density < 0.0);
    CHECK(!rep.passed);
}

TEST_CASE(butterfly_density_is_positive_for_healthy_slice) {
    const SviParams p{0.010, 0.30, -0.35, 0.00, 0.20};
    for (double k = -2.0; k <= 2.0; k += 0.05) {
        CHECK(butterfly_density(p, k) > 0.0);
    }
}

TEST_CASE(nelder_mead_minimises_rosenbrock) {
    const ObjectiveFn rosenbrock = [](const std::vector<double>& x) {
        const double a = 1.0 - x[0];
        const double b = x[1] - x[0] * x[0];
        return a * a + 100.0 * b * b;
    };
    NelderMeadOptions opt;
    opt.max_iterations = 10000;
    const NelderMeadResult r = nelder_mead(rosenbrock, {-1.2, 1.0}, {0.5, 0.5}, opt);
    CHECK(r.converged);
    CHECK_NEAR(r.x[0], 1.0, 0.0, 1e-4);
    CHECK_NEAR(r.x[1], 1.0, 0.0, 1e-4);
}
