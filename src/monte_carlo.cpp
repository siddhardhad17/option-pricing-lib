#include "optpricing/monte_carlo.hpp"

#include "optpricing/black_scholes.hpp"
#include "optpricing/rng.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <thread>
#include <utility>
#include <vector>

namespace optpricing {
namespace {

// Broadie, Glasserman & Kou (1997) continuity correction for discretely monitored barriers.
constexpr double kBgkShift = 0.5826;

void validate_config(const MonteCarloConfig& cfg, bool needs_steps) {
    if (cfg.paths < 4) {
        throw std::invalid_argument("paths must be >= 4");
    }
    if (needs_steps && cfg.steps == 0) {
        throw std::invalid_argument("steps must be >= 1");
    }
}

// Running sums of (y, x) pairs: y is the estimator, x is the control variate.
struct Accumulator {
    std::size_t n = 0;
    double sy = 0.0;
    double sx = 0.0;
    double syy = 0.0;
    double sxx = 0.0;
    double sxy = 0.0;

    void add(double y, double x) noexcept {
        ++n;
        sy += y;
        sx += x;
        syy += y * y;
        sxx += x * x;
        sxy += x * y;
    }

    void merge(const Accumulator& o) noexcept {
        n += o.n;
        sy += o.sy;
        sx += o.sx;
        syy += o.syy;
        sxx += o.sxx;
        sxy += o.sxy;
    }
};

unsigned resolve_threads(unsigned requested) {
    if (requested > 0) {
        return requested;
    }
    const unsigned hw = std::thread::hardware_concurrency();
    return hw == 0 ? 1u : hw;
}

std::size_t sample_count(const MonteCarloConfig& cfg) {
    return cfg.antithetic ? (cfg.paths + 1) / 2 : cfg.paths;
}

// Splits the samples across workers. Each worker owns a stream derived from (seed, worker id),
// so results depend only on the seed and the worker count.
template <class Sampler>
Accumulator run_parallel(std::size_t samples, unsigned requested_threads, std::uint64_t seed,
                         std::size_t buffer_size, const Sampler& sampler) {
    const std::size_t workers = std::min<std::size_t>(resolve_threads(requested_threads), samples);
    std::vector<Accumulator> parts(workers);

    auto work = [&](std::size_t id) {
        const std::size_t count = samples / workers + (id < samples % workers ? 1 : 0);
        NormalSource rng(derive_seed(seed, id));
        std::vector<double> buffer(buffer_size);
        Accumulator& acc = parts[id];
        for (std::size_t i = 0; i < count; ++i) {
            const std::pair<double, double> yx = sampler(rng, buffer);
            acc.add(yx.first, yx.second);
        }
    };

    if (workers == 1) {
        work(0);
    } else {
        std::vector<std::thread> pool;
        pool.reserve(workers);
        for (std::size_t id = 0; id < workers; ++id) {
            pool.emplace_back(work, id);
        }
        for (auto& t : pool) {
            t.join();
        }
    }

    Accumulator total;
    for (const auto& p : parts) {
        total.merge(p);
    }
    return total;
}

// Applies the control variate with the optimal coefficient beta = cov(y, x) / var(x),
// and returns the estimate with its standard error.
McEstimate finalize(const Accumulator& acc, double control_mean, bool use_control,
                    std::size_t total_paths) {
    const double n = static_cast<double>(acc.n);
    const double mean_y = acc.sy / n;
    const double mean_x = acc.sx / n;
    const double var_y = (acc.syy - n * mean_y * mean_y) / (n - 1.0);
    const double var_x = (acc.sxx - n * mean_x * mean_x) / (n - 1.0);
    const double cov_xy = (acc.sxy - n * mean_x * mean_y) / (n - 1.0);

    double beta = 0.0;
    if (use_control && var_x > 0.0) {
        beta = cov_xy / var_x;
    }
    const double price = mean_y - beta * (mean_x - control_mean);
    const double var_adj = std::max(0.0, var_y - 2.0 * beta * cov_xy + beta * beta * var_x);
    return {price, std::sqrt(var_adj / n), total_paths};
}

}  // namespace

McEstimate mc_european(const Market& m, const VanillaOption& o, double vol, const MonteCarloConfig& cfg) {
    validate_market(m);
    validate_vanilla(o);
    validate_vol(vol);
    validate_config(cfg, false);

    const double T = o.maturity;
    const double drift = (m.rate - m.dividend - 0.5 * vol * vol) * T;
    const double diffusion = vol * std::sqrt(T);
    const double discount = std::exp(-m.rate * T);
    const double phi = option_sign(o.type);
    const double strike = o.strike;
    const double spot = m.spot;
    const bool anti = cfg.antithetic;

    // Discounted payoff and discounted terminal price for one draw z.
    auto leg = [&](double z) -> std::pair<double, double> {
        const double st = spot * std::exp(drift + diffusion * z);
        return {discount * std::max(phi * (st - strike), 0.0), discount * st};
    };

    const auto sampler = [&](NormalSource& rng, std::vector<double>&) -> std::pair<double, double> {
        const double z = rng.next();
        const std::pair<double, double> first = leg(z);
        if (!anti) {
            return first;
        }
        const std::pair<double, double> second = leg(-z);
        return {0.5 * (first.first + second.first), 0.5 * (first.second + second.second)};
    };

    const std::size_t samples = sample_count(cfg);
    const Accumulator acc = run_parallel(samples, cfg.threads, cfg.seed, 1, sampler);
    // E[e^{-rT} S_T] = S e^{-qT} exactly, so the control has a known mean.
    const double control_mean = spot * std::exp(-m.dividend * T);
    return finalize(acc, control_mean, cfg.control_variate, samples * (anti ? 2 : 1));
}

McEstimate mc_barrier(const Market& m, const BarrierOption& o, double vol, const MonteCarloConfig& cfg) {
    validate_market(m);
    validate_barrier(o);
    validate_vol(vol);
    validate_config(cfg, true);

    const bool down = is_down_barrier(o.barrier);
    const bool knock_out = is_knock_out(o.barrier);
    const double T = o.maturity;
    const std::size_t steps = cfg.steps;
    const double dt = T / static_cast<double>(steps);
    const double drift_dt = (m.rate - m.dividend - 0.5 * vol * vol) * dt;
    const double diffusion_dt = vol * std::sqrt(dt);

    // Discrete monitoring misses some crossings, so the effective barrier moves toward the spot:
    // outward for down barriers (H rises), inward for up barriers (H falls).
    const double shift = kBgkShift * diffusion_dt;
    const double log_barrier = std::log(o.barrier_level) + (down ? shift : -shift);
    const bool breached_at_start = down ? (m.spot <= o.barrier_level) : (m.spot >= o.barrier_level);
    const double log_spot = std::log(m.spot);
    const double discount = std::exp(-m.rate * T);
    const double phi = option_sign(o.type);
    const double strike = o.strike;
    const bool anti = cfg.antithetic;

    // Returns (estimator, control) for the path driven by sign * z.
    auto path = [&](const std::vector<double>& z, double sign) -> std::pair<double, double> {
        double log_s = log_spot;
        bool touched = breached_at_start;
        for (std::size_t k = 0; k < steps; ++k) {
            log_s += drift_dt + diffusion_dt * sign * z[k];
            if (down ? (log_s <= log_barrier) : (log_s >= log_barrier)) {
                touched = true;
            }
        }
        const double vanilla_payoff = discount * std::max(phi * (std::exp(log_s) - strike), 0.0);
        const bool alive = knock_out ? !touched : touched;
        return {alive ? vanilla_payoff : 0.0, vanilla_payoff};
    };

    const auto sampler = [&](NormalSource& rng, std::vector<double>& z) -> std::pair<double, double> {
        for (double& v : z) {
            v = rng.next();
        }
        const std::pair<double, double> first = path(z, 1.0);
        if (!anti) {
            return first;
        }
        const std::pair<double, double> second = path(z, -1.0);
        return {0.5 * (first.first + second.first), 0.5 * (first.second + second.second)};
    };

    const std::size_t samples = sample_count(cfg);
    const Accumulator acc = run_parallel(samples, cfg.threads, cfg.seed, steps, sampler);
    // The control is the vanilla payoff on the same paths, whose mean is the exact BS price.
    const double control_mean = bs_price(m, VanillaOption{o.type, strike, T}, vol);
    return finalize(acc, control_mean, cfg.control_variate, samples * (anti ? 2 : 1));
}

}  // namespace optpricing
