#pragma once

namespace optpricing {

// Standard normal density.
double norm_pdf(double x) noexcept;

// Standard normal cumulative distribution, accurate in both tails (uses erfc).
double norm_cdf(double x) noexcept;

}  // namespace optpricing
