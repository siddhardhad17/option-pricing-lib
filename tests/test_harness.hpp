// Minimal self-registering test harness. No external dependencies.
#pragma once

#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace testing {

struct TestCase {
    std::string name;
    std::function<void()> fn;
};

inline std::vector<TestCase>& registry() {
    static std::vector<TestCase> tests;
    return tests;
}

inline int& failure_count() {
    static int failures = 0;
    return failures;
}

struct Registrar {
    Registrar(std::string name, std::function<void()> fn) {
        registry().push_back({std::move(name), std::move(fn)});
    }
};

inline void report(const char* file, int line, const std::string& what, const std::string& detail = "") {
    std::cerr << file << ":" << line << ": check failed: " << what;
    if (!detail.empty()) {
        std::cerr << " (" << detail << ")";
    }
    std::cerr << "\n";
    ++failure_count();
}

// |a - b| <= max(abs_tol, rel_tol * max(|a|, |b|))
inline bool near(double a, double b, double rel_tol, double abs_tol) {
    return std::fabs(a - b) <= std::max(abs_tol, rel_tol * std::max(std::fabs(a), std::fabs(b)));
}

}  // namespace testing

#define TEST_CASE(NAME)                                                       \
    static void NAME();                                                       \
    static const ::testing::Registrar NAME##_registrar(#NAME, &NAME);         \
    static void NAME()

#define CHECK(cond)                                                           \
    do {                                                                      \
        if (!(cond)) {                                                        \
            ::testing::report(__FILE__, __LINE__, #cond);                     \
        }                                                                     \
    } while (0)

#define CHECK_NEAR(actual, expected, rel_tol, abs_tol)                        \
    do {                                                                      \
        const double a_ = (actual);                                           \
        const double e_ = (expected);                                         \
        if (!::testing::near(a_, e_, (rel_tol), (abs_tol))) {                 \
            std::ostringstream os_;                                           \
            os_.precision(12);                                                \
            os_ << "actual=" << a_ << " expected=" << e_;                     \
            ::testing::report(__FILE__, __LINE__, #actual " ~ " #expected,    \
                              os_.str());                                     \
        }                                                                     \
    } while (0)

#define CHECK_THROWS_INVALID(expr)                                            \
    do {                                                                      \
        bool threw_ = false;                                                  \
        try {                                                                 \
            (void)(expr);                                                     \
        } catch (const std::invalid_argument&) {                              \
            threw_ = true;                                                    \
        }                                                                     \
        if (!threw_) {                                                        \
            ::testing::report(__FILE__, __LINE__,                             \
                              #expr " should throw std::invalid_argument");   \
        }                                                                     \
    } while (0)
