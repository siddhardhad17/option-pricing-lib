#include "optpricing/nelder_mead.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>

namespace optpricing {

NelderMeadResult nelder_mead(const ObjectiveFn& f, const std::vector<double>& x0,
                             const std::vector<double>& steps, const NelderMeadOptions& opt) {
    const std::size_t n = x0.size();
    if (n == 0 || steps.size() != n) {
        throw std::invalid_argument("nelder_mead: dimension mismatch");
    }

    using Point = std::vector<double>;
    constexpr double kReflect = 1.0;
    constexpr double kExpand = 2.0;
    constexpr double kContract = 0.5;
    constexpr double kShrink = 0.5;

    std::vector<Point> simplex(n + 1, x0);
    for (std::size_t i = 0; i < n; ++i) {
        simplex[i + 1][i] += steps[i];
    }
    std::vector<double> values(n + 1);
    for (std::size_t i = 0; i <= n; ++i) {
        values[i] = f(simplex[i]);
    }

    int iter = 0;
    bool converged = false;
    while (iter < opt.max_iterations) {
        ++iter;

        // Order vertices from best to worst.
        std::vector<std::size_t> order(n + 1);
        std::iota(order.begin(), order.end(), std::size_t{0});
        std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) { return values[a] < values[b]; });
        std::vector<Point> sorted_simplex(n + 1);
        std::vector<double> sorted_values(n + 1);
        for (std::size_t i = 0; i <= n; ++i) {
            sorted_simplex[i] = simplex[order[i]];
            sorted_values[i] = values[order[i]];
        }
        simplex.swap(sorted_simplex);
        values.swap(sorted_values);

        const double f_spread = values[n] - values[0];
        double x_spread = 0.0;
        for (std::size_t i = 1; i <= n; ++i) {
            for (std::size_t j = 0; j < n; ++j) {
                x_spread = std::max(x_spread, std::fabs(simplex[i][j] - simplex[0][j]));
            }
        }
        if (f_spread <= opt.f_tolerance * (1.0 + std::fabs(values[0])) && x_spread <= opt.x_tolerance) {
            converged = true;
            break;
        }

        // Centroid of every vertex except the worst.
        Point centroid(n, 0.0);
        for (std::size_t i = 0; i < n; ++i) {
            for (std::size_t j = 0; j < n; ++j) {
                centroid[j] += simplex[i][j];
            }
        }
        for (double& c : centroid) {
            c /= static_cast<double>(n);
        }

        const Point worst = simplex[n];
        // Point centroid + t * (centroid - worst).
        auto along = [&](double t) {
            Point p(n);
            for (std::size_t j = 0; j < n; ++j) {
                p[j] = centroid[j] + t * (centroid[j] - worst[j]);
            }
            return p;
        };

        const Point xr = along(kReflect);
        const double fr = f(xr);
        if (fr < values[0]) {
            const Point xe = along(kReflect * kExpand);
            const double fe = f(xe);
            if (fe < fr) {
                simplex[n] = xe;
                values[n] = fe;
            } else {
                simplex[n] = xr;
                values[n] = fr;
            }
        } else if (fr < values[n - 1]) {
            simplex[n] = xr;
            values[n] = fr;
        } else {
            bool accepted = false;
            if (fr < values[n]) {
                const Point xc = along(kReflect * kContract);
                const double fc = f(xc);
                if (fc <= fr) {
                    simplex[n] = xc;
                    values[n] = fc;
                    accepted = true;
                }
            } else {
                const Point xc = along(-kContract);
                const double fc = f(xc);
                if (fc < values[n]) {
                    simplex[n] = xc;
                    values[n] = fc;
                    accepted = true;
                }
            }
            if (!accepted) {
                for (std::size_t i = 1; i <= n; ++i) {
                    for (std::size_t j = 0; j < n; ++j) {
                        simplex[i][j] = simplex[0][j] + kShrink * (simplex[i][j] - simplex[0][j]);
                    }
                    values[i] = f(simplex[i]);
                }
            }
        }
    }

    std::size_t best = 0;
    for (std::size_t i = 1; i <= n; ++i) {
        if (values[i] < values[best]) {
            best = i;
        }
    }
    return {simplex[best], values[best], iter, converged};
}

}  // namespace optpricing
