#include "test_harness.hpp"

#include "optpricing/barrier.hpp"
#include "optpricing/black_scholes.hpp"

#include <vector>

using namespace optpricing;

namespace {

struct Case {
    BarrierType down_or_up;  // DownAndOut for down-side cases, UpAndOut for up-side cases
    double strike;
    double barrier;
};

}  // namespace

TEST_CASE(in_out_parity_holds_in_every_regime) {
    const Market m{100.0, 0.05, 0.02};
    const double vol = 0.30;
    const double T = 1.0;
    // Covers strike above and below the barrier, on both sides of spot.
    const std::vector<Case> cases = {
        {BarrierType::DownAndOut, 90.0, 85.0},   // K > H
        {BarrierType::DownAndOut, 80.0, 95.0},   // K <= H
        {BarrierType::DownAndOut, 90.0, 95.0},   // K <= H
        {BarrierType::DownAndOut, 100.0, 70.0},  // K > H
        {BarrierType::UpAndOut, 110.0, 115.0},   // K <= H
        {BarrierType::UpAndOut, 120.0, 105.0},   // K > H
        {BarrierType::UpAndOut, 100.0, 105.0},   // K <= H
    };
    for (const Case& c : cases) {
        for (OptionType type : {OptionType::Call, OptionType::Put}) {
            const VanillaOption v{type, c.strike, T};
            const bool down = c.down_or_up == BarrierType::DownAndOut;
            const BarrierOption out{type, c.down_or_up, c.strike, c.barrier, T};
            const BarrierOption in{type, down ? BarrierType::DownAndIn : BarrierType::UpAndIn, c.strike,
                                   c.barrier, T};
            const double sum = barrier_closed_form(m, out, vol) + barrier_closed_form(m, in, vol);
            CHECK_NEAR(sum, bs_price(m, v, vol), 0.0, 1e-9);
        }
    }
}

TEST_CASE(barrier_values_are_bounded_by_vanilla) {
    const Market m{100.0, 0.05, 0.0};
    const double vol = 0.25;
    for (OptionType type : {OptionType::Call, OptionType::Put}) {
        const double vanilla = bs_price(m, {type, 100.0, 1.0}, vol);
        const double do_price = barrier_closed_form(m, {type, BarrierType::DownAndOut, 100.0, 90.0, 1.0}, vol);
        const double uo_price = barrier_closed_form(m, {type, BarrierType::UpAndOut, 100.0, 115.0, 1.0}, vol);
        CHECK(do_price >= 0.0 && do_price <= vanilla + 1e-12);
        CHECK(uo_price >= 0.0 && uo_price <= vanilla + 1e-12);
    }
}

TEST_CASE(barrier_limits_recover_vanilla) {
    const Market m{100.0, 0.05, 0.02};
    const double vol = 0.3;
    const double T = 1.0;
    // A barrier that is practically unreachable must leave the vanilla price unchanged.
    for (OptionType type : {OptionType::Call, OptionType::Put}) {
        const double vanilla = bs_price(m, {type, 100.0, T}, vol);
        CHECK_NEAR(barrier_closed_form(m, {type, BarrierType::DownAndOut, 100.0, 1e-6, T}, vol), vanilla, 0.0, 1e-8);
        CHECK_NEAR(barrier_closed_form(m, {type, BarrierType::UpAndOut, 100.0, 1e8, T}, vol), vanilla, 0.0, 1e-8);
        CHECK_NEAR(barrier_closed_form(m, {type, BarrierType::DownAndIn, 100.0, 1e-6, T}, vol), 0.0, 0.0, 1e-8);
        CHECK_NEAR(barrier_closed_form(m, {type, BarrierType::UpAndIn, 100.0, 1e8, T}, vol), 0.0, 0.0, 1e-8);
    }
}

TEST_CASE(breached_barrier_at_inception) {
    const Market m{100.0, 0.05, 0.0};
    const double vol = 0.2;
    // Spot (100) is already at or below a down barrier at 105: knocked out at inception.
    const BarrierOption dout{OptionType::Call, BarrierType::DownAndOut, 95.0, 105.0, 1.0};
    const BarrierOption din{OptionType::Call, BarrierType::DownAndIn, 95.0, 105.0, 1.0};
    CHECK_NEAR(barrier_closed_form(m, dout, vol), 0.0, 0.0, 0.0);
    CHECK_NEAR(barrier_closed_form(m, din, vol), bs_price(m, {OptionType::Call, 95.0, 1.0}, vol), 0.0, 1e-12);
}

TEST_CASE(barrier_rejects_invalid_inputs) {
    const Market m{100.0, 0.05, 0.0};
    CHECK_THROWS_INVALID(barrier_closed_form(m, {OptionType::Call, BarrierType::DownAndOut, 100.0, 0.0, 1.0}, 0.2));
    CHECK_THROWS_INVALID(barrier_closed_form(m, {OptionType::Call, BarrierType::DownAndOut, 100.0, 90.0, 0.0}, 0.2));
    CHECK_THROWS_INVALID(barrier_closed_form(m, {OptionType::Call, BarrierType::DownAndOut, 100.0, 90.0, 1.0}, -1.0));
}
