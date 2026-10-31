#include "test_harness.hpp"

#include "optpricing/black_scholes.hpp"
#include "optpricing/greeks.hpp"

using namespace optpricing;

TEST_CASE(analytic_greeks_match_finite_differences) {
    const Market m{100.0, 0.04, 0.015};
    for (OptionType type : {OptionType::Call, OptionType::Put}) {
        for (double K : {80.0, 100.0, 120.0}) {
            for (double T : {0.5, 2.0}) {
                for (double vol : {0.15, 0.4}) {
                    const VanillaOption o{type, K, T};
                    const Greeks a = bs_greeks(m, o, vol);
                    const Greeks f = fd_greeks_vanilla(m, o, vol);
                    CHECK_NEAR(f.delta, a.delta, 1e-6, 1e-7);
                    CHECK_NEAR(f.gamma, a.gamma, 1e-4, 1e-6);
                    CHECK_NEAR(f.vega, a.vega, 1e-6, 1e-6);
                    CHECK_NEAR(f.theta, a.theta, 1e-6, 1e-6);
                    CHECK_NEAR(f.rho, a.rho, 1e-6, 1e-6);
                }
            }
        }
    }
}

TEST_CASE(greek_signs_are_sensible) {
    const Market m{100.0, 0.03, 0.0};
    const VanillaOption call{OptionType::Call, 100.0, 1.0};
    const VanillaOption put{OptionType::Put, 100.0, 1.0};
    const Greeks gc = bs_greeks(m, call, 0.2);
    const Greeks gp = bs_greeks(m, put, 0.2);
    CHECK(gc.delta > 0.0 && gc.delta < 1.0);
    CHECK(gp.delta < 0.0 && gp.delta > -1.0);
    CHECK(gc.gamma > 0.0 && gp.gamma > 0.0);
    CHECK(gc.vega > 0.0 && gp.vega > 0.0);
    CHECK(gc.rho > 0.0 && gp.rho < 0.0);
    CHECK(gc.theta < 0.0);  // long at-the-money calls lose time value
    CHECK_NEAR(gc.delta - gp.delta, 1.0, 1e-12, 1e-12);  // put-call delta parity with zero dividend
    CHECK_NEAR(gc.gamma, gp.gamma, 1e-12, 1e-14);
    CHECK_NEAR(gc.vega, gp.vega, 1e-12, 1e-12);
}
