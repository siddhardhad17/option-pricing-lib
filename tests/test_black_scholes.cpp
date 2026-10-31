#include "test_harness.hpp"

#include "optpricing/black_scholes.hpp"

#include <cmath>

using namespace optpricing;

TEST_CASE(bs_matches_reference_values) {
    // Hull, "Options, Futures and Other Derivatives": S = 42, K = 40, r = 10%, sigma = 20%, T = 0.5.
    const Market m{42.0, 0.10, 0.0};
    CHECK_NEAR(bs_price(m, {OptionType::Call, 40.0, 0.5}, 0.20), 4.7594, 0.0, 1e-4);
    CHECK_NEAR(bs_price(m, {OptionType::Put, 40.0, 0.5}, 0.20), 0.8086, 0.0, 1e-4);
}

TEST_CASE(put_call_parity_holds_with_dividends) {
    const Market m{100.0, 0.05, 0.03};
    for (double K : {70.0, 100.0, 130.0}) {
        for (double vol : {0.05, 0.25, 0.8}) {
            const double T = 1.3;
            const double call = bs_price(m, {OptionType::Call, K, T}, vol);
            const double put = bs_price(m, {OptionType::Put, K, T}, vol);
            const double parity = m.spot * std::exp(-m.dividend * T) - K * std::exp(-m.rate * T);
            CHECK_NEAR(call - put, parity, 0.0, 1e-10);
        }
    }
}

TEST_CASE(implied_vol_round_trip) {
    // Vols below ~10% combined with deep in-the-money strikes give time values under 1e-12,
    // which double precision cannot resolve, so those volatilities are not identifiable.
    const Market m{100.0, 0.02, 0.01};
    for (double vol : {0.1, 0.2, 0.6, 1.5}) {
        for (double K : {70.0, 100.0, 140.0}) {
            for (OptionType type : {OptionType::Call, OptionType::Put}) {
                const VanillaOption o{type, K, 0.75};
                const double price = bs_price(m, o, vol);
                CHECK_NEAR(bs_implied_vol(m, o, price), vol, 1e-7, 1e-8);
            }
        }
    }
}

TEST_CASE(implied_vol_rejects_arbitrage_prices) {
    const Market m{100.0, 0.05, 0.0};
    const VanillaOption call{OptionType::Call, 90.0, 1.0};
    // Below intrinsic value (discounted) and above the forward bound must both be rejected.
    CHECK_THROWS_INVALID(bs_implied_vol(m, call, 0.0));
    CHECK_THROWS_INVALID(bs_implied_vol(m, call, 150.0));
}

TEST_CASE(invalid_inputs_are_rejected) {
    const Market m{100.0, 0.05, 0.0};
    const VanillaOption o{OptionType::Call, 100.0, 1.0};
    CHECK_THROWS_INVALID(bs_price(m, o, 0.0));
    CHECK_THROWS_INVALID(bs_price(m, o, -0.2));
    CHECK_THROWS_INVALID(bs_price(Market{-1.0, 0.0, 0.0}, o, 0.2));
    CHECK_THROWS_INVALID(bs_price(m, VanillaOption{OptionType::Call, 100.0, 0.0}, 0.2));
    CHECK_THROWS_INVALID(bs_price(m, VanillaOption{OptionType::Call, std::nan(""), 1.0}, 0.2));
}
