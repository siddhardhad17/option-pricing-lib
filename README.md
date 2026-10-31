# optpricing: option pricing and SVI volatility surface library

A C++17 library for European and barrier option pricing, Greeks, Monte Carlo variance reduction,
SVI volatility-surface calibration, and static no-arbitrage checks. It has no third-party
dependencies, builds with CMake, and is tested against closed-form benchmarks.

## Features

| Area | What is implemented |
|---|---|
| Vanilla pricing | Black-Scholes-Merton with continuous dividend yield; implied volatility by bisection |
| Greeks | Closed-form delta, gamma, vega, theta, rho; generic central finite differences for any pricer |
| Barrier pricing | Closed-form continuously monitored barriers (down/up, in/out, call/put, zero rebate) |
| Monte Carlo | GBM paths for European options; discretely monitored barriers with the Broadie-Glasserman-Kou correction |
| Variance reduction | Antithetic variates; control variates (discounted terminal price for European, vanilla payoff for barrier) |
| Parallelism | `std::thread` workers, each with a reproducible stream derived from `(seed, worker id)` |
| Volatility surface | Raw SVI slice fitting (Nelder-Mead with constraints and restarts) from prices or implied vols |
| Static arbitrage | Butterfly density g(k) >= 0, Roger Lee wing bound, calendar monotonicity, positive variance |

## Repository layout

```
include/optpricing/   public headers (the API)
src/                  implementation
tests/                dependency-free unit tests (27 cases)
apps/                 price_demo, benchmark, svi_demo
cmake/                package config template for find_package
.github/workflows/    CI: GCC and Clang builds, sanitizer build
```

## Build and test

Requirements: a C++17 compiler (tested with GCC 13 and Clang 18 on Linux), CMake 3.16+, and a threads
library. It has not been tested on MSVC or macOS.

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

Useful options:

| Option | Default | Effect |
|---|---|---|
| `OPTPRICING_BUILD_TESTS` | `ON` | Build the unit tests |
| `OPTPRICING_BUILD_APPS` | `ON` | Build `price_demo`, `benchmark`, `svi_demo` |
| `OPTPRICING_SANITIZE` | `OFF` | Build with AddressSanitizer and UndefinedBehaviorSanitizer |
| `OPTPRICING_WARNINGS_AS_ERRORS` | `OFF` | Add `-Werror` (CI enables it) |

## Install and use from another CMake project

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/opt/optpricing
cmake --build build && cmake --install build
```

```cmake
find_package(optpricing 1.0 REQUIRED)
target_link_libraries(my_app PRIVATE optpricing::optpricing)
```

Example program:

```cpp
#include <optpricing/barrier.hpp>
#include <optpricing/black_scholes.hpp>
#include <optpricing/greeks.hpp>
#include <optpricing/monte_carlo.hpp>

#include <iostream>

using namespace optpricing;

int main() {
    const Market m{100.0, 0.05, 0.02};  // spot, risk-free rate, dividend yield (continuous)
    const double vol = 0.25;

    const VanillaOption call{OptionType::Call, 100.0, 1.0};
    std::cout << "BS call price : " << bs_price(m, call, vol) << "\n";
    std::cout << "BS call delta : " << bs_greeks(m, call, vol).delta << "\n";

    const BarrierOption dout{OptionType::Call, BarrierType::DownAndOut, 100.0, 85.0, 1.0};
    std::cout << "DO call closed form : " << barrier_closed_form(m, dout, vol) << "\n";

    MonteCarloConfig cfg;   // defaults: 100k paths, 252 steps, antithetic + control variate
    cfg.paths = 1000000;
    cfg.threads = 0;        // 0 = hardware concurrency
    const McEstimate mc = mc_barrier(m, dout, vol, cfg);
    std::cout << "DO call Monte Carlo : " << mc.price << " +/- " << mc.std_error << "\n";
}
```

## Conventions

- Rates are continuously compounded and annualised; maturities are in years.
- `Greeks::vega` is per 1.00 of volatility (divide by 100 for a 1-vol-point move).
- `Greeks::theta` is `dV/dt = -dV/dT`, per year.
- `Greeks::rho` is per 1.00 of the rate.
- Invalid inputs (non-positive spot, strike, maturity or volatility, non-finite numbers) throw
  `std::invalid_argument`.
- Implied volatility raises `std::invalid_argument` when the price is at or outside the
  no-arbitrage bounds. A price within rounding error of intrinsic value carries no identifiable
  volatility, so it is rejected rather than guessed.

## Method notes

**Closed-form barriers.** The formulas are the Reiner-Rubinstein terms A, B, C, D with the case table
for each barrier type, split by whether the strike lies above the barrier. The case table was checked
cell by cell against a 400k-path Monte Carlo before it was coded, and every knock-in/knock-out pair
satisfies in-out parity. Knock-outs are zero when the barrier is breached at inception.

**Monte Carlo.** European options use the exact terminal draw under GBM. Barrier options use
`steps` monitoring dates. The barrier is moved toward the spot by the BGK shift
`exp(±0.5826 sigma sqrt(dt))`, which removes most of the discrete-monitoring bias. Standard errors come
from the sample variance of the estimator, so they stay valid with antithetic pairs and control variates.

**Reproducibility.** The generator is Box-Muller over `std::mt19937_64`, not `std::normal_distribution`,
so the stream does not depend on the standard library. A fixed `(seed, threads, paths, antithetic)` gives
bit-identical results. Changing the thread count changes the stream partitioning and therefore the result
within statistical error. With `threads = 0` the result also depends on the machine's core count.

**SVI.** The raw parameterisation `w(k) = a + b(rho(k - m) + sqrt((k - m)^2 + sigma^2))` is fitted to
total variance `sigma_imp^2 T` against `k = ln(K / F)`. Constraints (`b >= 0`, `|rho| < 1`, `sigma > 0`,
non-negative minimum variance) are enforced by penalty. Nelder-Mead runs from a data-driven start, then
restarts from the incumbent optimum.

**Static checks.** `check_static_arbitrage` reports:

- *Butterfly:* Durrleman's density `g(k) >= 0` on a grid. This is the condition for a non-negative
  risk-neutral density.
- *Calendar:* `w(T2, k) >= w(T1, k)` for consecutive expiries on the grid.
- *Lee bound:* `b(1 +/- rho) <= 2`, the necessary condition for finite moments.
- *Positivity:* `w(k) > 0` on the grid.

These are necessary conditions evaluated on a finite grid. They do not prove the absence of
arbitrage everywhere, and they do not check put-call consistency of raw quotes. Quotes should be
cleaned and checked upstream.

## Verification

Tests are in `tests/` and run with `ctest`. They cover:

- Black-Scholes against the Hull reference values, put-call parity with dividends, and implied-volatility round trips.
- Analytic Greeks against central finite differences for calls and puts over several strikes and maturities.
- In-out parity for all eight barrier cells in both strike regimes, the vanilla limits, and breached barriers.
- Monte Carlo against closed forms (European and eight barrier cells), variance-reduction effects, seed reproducibility, and thread-count consistency.
- SVI identities, recovery of a synthetic slice from prices, Nelder-Mead on Rosenbrock, and detection of butterfly, calendar and Lee-bound violations.

Also checked: clean builds with GCC and Clang under `-Wall -Wextra -Wpedantic -Wshadow -Werror`;
AddressSanitizer and UndefinedBehaviorSanitizer runs with no reports; ThreadSanitizer runs with no
data-race reports; and a test consumer project that builds against the installed package.

Sample `benchmark` output (2-core VM, S = 100, K = 100, sigma = 0.25, T = 1, r = 5%, q = 2%):

| Option | Paths | Variant | Abs. error | Std. error | Std. error reduction |
|---|---|---|---|---|---|
| call | 1,000,000 | plain | 2.9e-03 | 1.75e-02 | 1.00x |
| call | 1,000,000 | antithetic | 3.2e-03 | 1.36e-02 | 1.29x |
| call | 1,000,000 | control variate | 1.2e-03 | 7.10e-03 | 2.47x |
| call | 1,000,000 | both | 9.0e-04 | 3.55e-03 | 4.95x |

Barrier Monte Carlo (252 steps, 200k paths, both variance reductions) against the closed form: all eight
cells have |z| < 1 (for example down-and-out call 9.9153 vs 9.9234, standard error 0.011).

Scaling on the same VM (2 hardware threads). Wall times vary by about 10-20% between runs on shared
hardware:

- European, one thread: about 56 million paths per second, flat across 0.25M to 4M paths.
- European, 4M paths: 1 thread 0.10 s, 2 threads 0.08 s, 4 threads 0.05 s, 8 threads 0.04 s.
- Barrier, 400k paths, 252 steps: 2 threads give close to 2x. Beyond 2 threads there is no
  further gain on this machine, as expected with only 2 cores.

Run `build/benchmark out_dir` on your own hardware for numbers that apply to your machine.

## Limitations

- Barrier closed forms cover continuous monitoring with zero rebate. Discretely monitored barriers,
  rebates, and double barriers are not implemented.
- Only European and barrier options are priced. There is no American exercise, local volatility,
  or Heston model.
- SVI is fitted to one expiry at a time. The surface-level calendar check runs on the fitted slices.
- Static checks are necessary conditions on a grid, not a proof of arbitrage-freeness.
- Monte Carlo barrier prices carry a small residual bias from finite monitoring steps (roughly 0.1-0.2%
  of price at 252 steps in the runs above). Increase `steps` or use a Brownian-bridge correction for tighter accuracy.

## Continuous integration

`.github/workflows/ci.yml` builds with GCC and Clang using `-Werror`, runs the test suite, and runs a
separate Debug build with AddressSanitizer and UBSan.

## License

MIT. See [LICENSE](LICENSE).
