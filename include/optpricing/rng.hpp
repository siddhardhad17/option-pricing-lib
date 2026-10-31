#pragma once

#include <cmath>
#include <cstdint>
#include <random>

namespace optpricing {

// SplitMix64 mixing step, used to derive independent per-thread seeds.
inline std::uint64_t splitmix64(std::uint64_t x) noexcept {
    x += 0x9E3779B97F4A7C15ULL;
    x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
    x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
    return x ^ (x >> 31);
}

inline std::uint64_t derive_seed(std::uint64_t seed, std::uint64_t stream) noexcept {
    return splitmix64(seed ^ splitmix64(stream + 1));
}

// Standard normal generator (Box-Muller over mt19937_64). Unlike std::normal_distribution,
// the output stream is identical across standard library implementations.
class NormalSource {
public:
    explicit NormalSource(std::uint64_t seed) : engine_(seed) {}

    double next() noexcept {
        if (has_spare_) {
            has_spare_ = false;
            return spare_;
        }
        const double u1 = uniform_open();
        const double u2 = uniform_open();
        const double r = std::sqrt(-2.0 * std::log(u1));
        const double theta = 2.0 * 3.14159265358979323846 * u2;
        spare_ = r * std::sin(theta);
        has_spare_ = true;
        return r * std::cos(theta);
    }

private:
    // Uniform on the open interval (0, 1) with 53 bits of resolution.
    double uniform_open() noexcept {
        return (static_cast<double>(engine_() >> 11) + 0.5) * 0x1.0p-53;
    }

    std::mt19937_64 engine_;
    double spare_ = 0.0;
    bool has_spare_ = false;
};

}  // namespace optpricing
