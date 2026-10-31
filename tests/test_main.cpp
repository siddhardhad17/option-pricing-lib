#include "test_harness.hpp"

#include <exception>

int main() {
    int failed_tests = 0;
    for (const auto& test : testing::registry()) {
        const int before = testing::failure_count();
        try {
            test.fn();
        } catch (const std::exception& e) {
            testing::report(__FILE__, __LINE__, test.name, std::string("unexpected exception: ") + e.what());
        }
        const bool ok = testing::failure_count() == before;
        std::cout << (ok ? "[  PASS  ] " : "[  FAIL  ] ") << test.name << "\n";
        if (!ok) {
            ++failed_tests;
        }
    }
    std::cout << "\n" << (testing::registry().size() - static_cast<std::size_t>(failed_tests)) << "/"
              << testing::registry().size() << " tests passed\n";
    return failed_tests == 0 ? 0 : 1;
}
