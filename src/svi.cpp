#include "optpricing/svi.hpp"

#include "optpricing/black_scholes.hpp"
#include "optpricing/nelder_mead.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace optpricing {
namespace {

// Amount by which parameters violate the SVI constraints (0 when valid).
double constraint_violation(const SviParams& p) {
    if (!std::isfinite(p.a) || !std::isfinite(p.b) || !std::isfinite(p.rho) || !std::isfinite(p.m) ||
        !std::isfinite(p.sigma)) {
        return std::numeric_limits<double>::infinity();
    }
    double v = 0.0;
    v += std::max(0.0, -p.b);
    v += std::max(0.0, -p.sigma);
    v += std::max(0.0, std::fabs(p.rho) - 0.999);
    const double rho_c = std::clamp(p.rho, -0.999, 0.999);
    v += std::max(0.0, -(p.a + p.b * p.sigma * std::sqrt(1.0 - rho_c * rho_c)));
    return v;
}

}  // namespace

double svi_total_variance(const SviParams& p, double k) noexcept {
    const double x = k - p.m;
    return p.a + p.b * (p.rho * x + std::sqrt(x * x + p.sigma * p.sigma));
}

double svi_total_variance_slope(const SviParams& p, double k) noexcept {
    const double x = k - p.m;
    return p.b * (p.rho + x / std::sqrt(x * x + p.sigma * p.sigma));
}

double svi_total_variance_curvature(const SviParams& p, double k) noexcept {
    const double x = k - p.m;
    const double r2 = x * x + p.sigma * p.sigma;
    return p.b * p.sigma * p.sigma / (r2 * std::sqrt(r2));
}

bool svi_params_valid(const SviParams& p) noexcept {
    return constraint_violation(p) == 0.0 && p.sigma > 0.0 && std::fabs(p.rho) < 1.0;
}

SmileQuote quote_from_price(const Market& m, const VanillaOption& o, double price) {
    return {o.strike, bs_implied_vol(m, o, price)};
}

SviFit fit_svi_slice(const std::vector<SmileQuote>& quotes, double forward, double maturity,
                     const SviFitOptions& opt) {
    if (quotes.size() < 5) {
        throw std::invalid_argument("fit_svi_slice needs at least five quotes");
    }
    if (!(forward > 0.0) || !(maturity > 0.0)) {
        throw std::invalid_argument("forward and maturity must be > 0");
    }

    const std::size_t n = quotes.size();
    std::vector<double> k(n);
    std::vector<double> w(n);
    std::vector<double> iv(n);
    for (std::size_t i = 0; i < n; ++i) {
        if (!(quotes[i].strike > 0.0) || !(quotes[i].implied_vol > 0.0)) {
            throw std::invalid_argument("quotes need positive strikes and implied volatilities");
        }
        k[i] = std::log(quotes[i].strike / forward);
        iv[i] = quotes[i].implied_vol;
        w[i] = iv[i] * iv[i] * maturity;
    }

    // Initial guess: vertex at the observed minimum of total variance.
    std::size_t imin = 0;
    for (std::size_t i = 1; i < n; ++i) {
        if (w[i] < w[imin]) {
            imin = i;
        }
    }
    const std::vector<double> x0 = {0.5 * w[imin], 0.1, -0.3, k[imin], 0.1};
    const std::vector<double> first_steps = {0.02, 0.05, 0.1, 0.05, 0.05};
    const std::vector<double> restart_steps = {0.005, 0.01, 0.02, 0.01, 0.01};

    const ObjectiveFn objective = [&](const std::vector<double>& x) {
        const SviParams p{x[0], x[1], x[2], x[3], x[4]};
        const double violation = constraint_violation(p);
        if (violation > 0.0) {
            return 1e3 * (1.0 + violation);
        }
        double sse = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            const double d = svi_total_variance(p, k[i]) - w[i];
            sse += d * d;
        }
        return sse / static_cast<double>(n);
    };

    NelderMeadOptions nm;
    nm.max_iterations = opt.max_iterations;
    nm.f_tolerance = opt.tolerance;
    nm.x_tolerance = 1e-12;

    NelderMeadResult best = nelder_mead(objective, x0, first_steps, nm);
    int iterations = best.iterations;
    for (int r = 0; r < opt.restarts; ++r) {
        NelderMeadResult trial = nelder_mead(objective, best.x, restart_steps, nm);
        iterations += trial.iterations;
        if (!(trial.fx < best.fx)) {
            break;
        }
        best = trial;
    }

    SviFit fit{};
    fit.params = SviParams{best.x[0], best.x[1], best.x[2], best.x[3], best.x[4]};
    fit.iterations = iterations;
    fit.converged = best.converged && constraint_violation(fit.params) == 0.0;

    double se_w = 0.0;
    double se_iv = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        const double model_w = svi_total_variance(fit.params, k[i]);
        se_w += (model_w - w[i]) * (model_w - w[i]);
        const double model_iv = std::sqrt(std::max(model_w, 0.0) / maturity);
        se_iv += (model_iv - iv[i]) * (model_iv - iv[i]);
    }
    fit.rmse_total_variance = std::sqrt(se_w / static_cast<double>(n));
    fit.rmse_implied_vol = std::sqrt(se_iv / static_cast<double>(n));
    return fit;
}

}  // namespace optpricing
