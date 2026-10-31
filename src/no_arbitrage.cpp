#include "optpricing/no_arbitrage.hpp"

#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>

namespace optpricing {
namespace {

constexpr double kTolerance = 1e-12;

std::vector<double> linspace(double lo, double hi, int points) {
    std::vector<double> out(static_cast<std::size_t>(points));
    const double step = (hi - lo) / static_cast<double>(points - 1);
    for (int i = 0; i < points; ++i) {
        out[static_cast<std::size_t>(i)] = lo + step * static_cast<double>(i);
    }
    return out;
}

std::string format(double x) {
    std::ostringstream os;
    os.precision(6);
    os << x;
    return os.str();
}

}  // namespace

double butterfly_density(const SviParams& p, double k) noexcept {
    const double w = svi_total_variance(p, k);
    const double w1 = svi_total_variance_slope(p, k);
    const double w2 = svi_total_variance_curvature(p, k);
    const double t = 1.0 - k * w1 / (2.0 * w);
    return t * t - 0.25 * w1 * w1 * (1.0 / w + 0.25) + 0.5 * w2;
}

ArbitrageReport check_static_arbitrage(const std::vector<SviSlice>& slices, const CheckGrid& grid) {
    if (slices.empty()) {
        throw std::invalid_argument("check_static_arbitrage needs at least one slice");
    }
    if (grid.points < 2 || !(grid.k_max > grid.k_min)) {
        throw std::invalid_argument("invalid log-moneyness grid");
    }
    for (std::size_t i = 1; i < slices.size(); ++i) {
        if (!(slices[i].maturity > slices[i - 1].maturity)) {
            throw std::invalid_argument("slices must have strictly increasing maturities");
        }
    }

    ArbitrageReport rep;
    const std::vector<double> ks = linspace(grid.k_min, grid.k_max, grid.points);

    for (std::size_t s = 0; s < slices.size(); ++s) {
        const SviParams& p = slices[s].params;
        const std::string label = "slice T=" + format(slices[s].maturity);

        if (!svi_params_valid(p)) {
            rep.passed = rep.parameters_valid = false;
            rep.violations.push_back(label + ": parameters violate SVI constraints");
            continue;
        }

        // Roger Lee moment bound: asymptotic wing slopes of total variance must not exceed 2.
        if (p.b * (1.0 + p.rho) > 2.0 + kTolerance || p.b * (1.0 - p.rho) > 2.0 + kTolerance) {
            rep.passed = rep.lee_wing_bounds = false;
            rep.violations.push_back(label + ": wing slope b(1 +/- rho) exceeds Lee bound 2");
        }

        double min_g = std::numeric_limits<double>::infinity();
        double min_g_k = ks.front();
        bool positive = true;
        for (double k : ks) {
            if (!(svi_total_variance(p, k) > 0.0)) {
                positive = false;
                continue;
            }
            const double g = butterfly_density(p, k);
            if (g < min_g) {
                min_g = g;
                min_g_k = k;
            }
        }
        if (!positive) {
            rep.passed = rep.positive_variance = false;
            rep.violations.push_back(label + ": total variance is non-positive on the grid");
        }
        if (std::isfinite(min_g)) {
            rep.min_butterfly_density = std::min(rep.min_butterfly_density, min_g);
            if (min_g < -kTolerance) {
                rep.passed = rep.butterfly = false;
                rep.violations.push_back(label + ": butterfly density negative (g = " + format(min_g) +
                                         ") at k = " + format(min_g_k));
            }
        }
    }

    // Calendar spread: total variance must not decrease with maturity at any fixed log-moneyness.
    for (std::size_t s = 1; s < slices.size(); ++s) {
        const SviParams& shorter = slices[s - 1].params;
        const SviParams& longer = slices[s].params;
        if (!svi_params_valid(shorter) || !svi_params_valid(longer)) {
            continue;
        }
        double worst = std::numeric_limits<double>::infinity();
        double worst_k = ks.front();
        for (double k : ks) {
            const double gap = svi_total_variance(longer, k) - svi_total_variance(shorter, k);
            if (gap < worst) {
                worst = gap;
                worst_k = k;
            }
        }
        rep.worst_calendar_gap = std::min(rep.worst_calendar_gap, worst);
        if (worst < -kTolerance) {
            rep.passed = rep.calendar = false;
            rep.violations.push_back("calendar arbitrage between T=" + format(slices[s - 1].maturity) + " and T=" +
                                     format(slices[s].maturity) + " at k = " + format(worst_k));
        }
    }
    return rep;
}

}  // namespace optpricing
