#include "test_harness.hpp"

#include "optpricing/barrier.hpp"
#include "optpricing/black_scholes.hpp"
#include "optpricing/monte_carlo.hpp"

#include <cmath>

using namespace optpricing;

namespace {

MonteCarloConfig config(std::size_t paths, bool anti, bool cv, unsigned threads = 1) {
    MonteCarloConfig cfg;
    cfg.paths = paths;
    cfg.antithetic = anti;
    cfg.control_variate = cv;
    cfg.threads = threads;
    cfg.seed = 7;
    return cfg;
}

}  // namespace

TEST_CASE(european_mc_agrees_with_black_scholes) {
    const Market m{100.0, 0.05, 0.01};
    const double vol = 0.25;
    for (OptionType type : {OptionType::Call, OptionType::Put}) {
        const VanillaOption o{type, 105.0, 1.0};
        const double exact = bs_price(m, o, vol);
        const McEstimate mc = mc_european(m, o, vol, config(400000, true, true));
        CHECK(std::fabs(mc.price - exact) <= 4.0 * mc.std_error);
        CHECK(mc.std_error < 0.01);
        CHECK(mc.paths == 400000);
    }
}

TEST_CASE(variance_reduction_lowers_standard_error) {
    const Market m{100.0, 0.05, 0.0};
    const double vol = 0.2;
    const VanillaOption call{OptionType::Call, 100.0, 1.0};
    const double plain = mc_european(m, call, vol, config(200000, false, false)).std_error;
    const double anti = mc_european(m, call, vol, config(200000, true, false)).std_error;
    const double cv = mc_european(m, call, vol, config(200000, false, true)).std_error;
    const double both = mc_european(m, call, vol, config(200000, true, true)).std_error;
    CHECK(anti < plain);
    CHECK(cv < 0.5 * plain);
    CHECK(both <= cv);
}

TEST_CASE(fixed_seed_is_reproducible) {
    const Market m{100.0, 0.05, 0.0};
    const VanillaOption o{OptionType::Put, 95.0, 0.5};
    const McEstimate a = mc_european(m, o, 0.3, config(50000, true, true, 2));
    const McEstimate b = mc_european(m, o, 0.3, config(50000, true, true, 2));
    CHECK(a.price == b.price);
    CHECK(a.std_error == b.std_error);
}

TEST_CASE(multithreaded_estimate_agrees_with_single_thread) {
    const Market m{100.0, 0.05, 0.0};
    const VanillaOption o{OptionType::Call, 110.0, 1.0};
    const double vol = 0.3;
    const double exact = bs_price(m, o, vol);
    const McEstimate one = mc_european(m, o, vol, config(300000, true, true, 1));
    const McEstimate four = mc_european(m, o, vol, config(300000, true, true, 4));
    CHECK(std::fabs(four.price - exact) <= 4.0 * four.std_error);
    CHECK(std::fabs(one.price - four.price) <= 4.0 * std::sqrt(one.std_error * one.std_error +
                                                              four.std_error * four.std_error));
}

TEST_CASE(barrier_mc_agrees_with_closed_form) {
    const Market m{100.0, 0.05, 0.02};
    const double vol = 0.3;
    const double T = 1.0;
    struct Cell {
        BarrierType type;
        double strike;
        double barrier;
    };
    const Cell cells[] = {
        {BarrierType::DownAndOut, 100.0, 85.0},
        {BarrierType::DownAndIn, 90.0, 80.0},
        {BarrierType::UpAndOut, 100.0, 118.0},
        {BarrierType::UpAndIn, 110.0, 120.0},
    };
    for (const Cell& c : cells) {
        for (OptionType type : {OptionType::Call, OptionType::Put}) {
            const BarrierOption o{type, c.type, c.strike, c.barrier, T};
            const double exact = barrier_closed_form(m, o, vol);
            MonteCarloConfig cfg = config(200000, true, true);
            cfg.steps = 252;
            const McEstimate mc = mc_barrier(m, o, vol, cfg);
            // Residual discretisation bias after the BGK correction is well below 1% of price here.
            CHECK(std::fabs(mc.price - exact) <= 4.0 * mc.std_error + 0.01 * exact + 1e-3);
        }
    }
}

TEST_CASE(barrier_mc_rejects_zero_steps) {
    const Market m{100.0, 0.05, 0.0};
    MonteCarloConfig cfg = config(1000, true, true);
    cfg.steps = 0;
    CHECK_THROWS_INVALID(mc_barrier(m, {OptionType::Call, BarrierType::DownAndOut, 100.0, 90.0, 1.0}, 0.2, cfg));
}
