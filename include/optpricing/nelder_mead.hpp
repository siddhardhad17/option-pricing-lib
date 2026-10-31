#pragma once

#include <cstddef>
#include <functional>
#include <vector>

namespace optpricing {

using ObjectiveFn = std::function<double(const std::vector<double>&)>;

struct NelderMeadOptions {
    int max_iterations = 20000;
    double f_tolerance = 1e-12;  // spread of objective values across the simplex
    double x_tolerance = 1e-10;  // largest vertex distance from the best vertex
};

struct NelderMeadResult {
    std::vector<double> x;  // best point found
    double fx;              // objective at x
    int iterations;
    bool converged;
};

// Derivative-free simplex minimiser (Nelder & Mead, 1965) with standard coefficients.
// `steps[i]` is the initial simplex offset along coordinate i.
NelderMeadResult nelder_mead(const ObjectiveFn& f, const std::vector<double>& x0,
                             const std::vector<double>& steps, const NelderMeadOptions& opt = {});

}  // namespace optpricing
