// Prices one European and one barrier option, compares closed form, Monte Carlo and Greeks.
//
// Usage: price_demo [--spot S] [--strike K] [--barrier H] [--rate r] [--div q] [--vol s]
//                   [--maturity T] [--paths N] [--steps M] [--threads P] [--seed X]

#include "optpricing/barrier.hpp"
#include "optpricing/black_scholes.hpp"
#include "optpricing/greeks.hpp"
#include "optpricing/monte_carlo.hpp"

#include <cstdlib>
#include <exception>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>

using namespace optpricing;

namespace {

void print_usage() {
    std::cout << "Usage: price_demo [--spot S] [--strike K] [--barrier H] [--rate r] [--div q]\n"
                 "                  [--vol s] [--maturity T] [--paths N] [--steps M]\n"
                 "                  [--threads P] [--seed X]\n";
}

double as_double(const std::string& key, const std::string& value) {
    try {
        std::size_t used = 0;
        const double v = std::stod(value, &used);
        if (used != value.size()) {
            throw std::invalid_argument("trailing characters");
        }
        return v;
    } catch (const std::exception&) {
        throw std::invalid_argument("invalid number for " + key + ": " + value);
    }
}

}  // namespace

int main(int argc, char** argv) {
    std::map<std::string, std::string> args;
    try {
        for (int i = 1; i < argc; ++i) {
            const std::string key = argv[i];
            if (key == "--help" || key == "-h") {
                print_usage();
                return 0;
            }
            if (key.rfind("--", 0) != 0 || i + 1 >= argc) {
                throw std::invalid_argument("expected --name value pairs, got: " + key);
            }
            args[key.substr(2)] = argv[++i];
        }
        auto get = [&](const char* name, double fallback) {
            const auto it = args.find(name);
            return it == args.end() ? fallback : as_double(name, it->second);
        };

        const Market market{get("spot", 100.0), get("rate", 0.05), get("div", 0.02)};
        const double strike = get("strike", 100.0);
        const double barrier = get("barrier", 85.0);
        const double vol = get("vol", 0.25);
        const double maturity = get("maturity", 1.0);

        MonteCarloConfig cfg;
        cfg.paths = static_cast<std::size_t>(get("paths", 400000));
        cfg.steps = static_cast<std::size_t>(get("steps", 252));
        cfg.threads = static_cast<unsigned>(get("threads", 0));
        cfg.seed = static_cast<std::uint64_t>(get("seed", 20260101));

        std::cout << std::fixed << std::setprecision(6);
        std::cout << "Market: S=" << market.spot << " r=" << market.rate << " q=" << market.dividend
                  << " sigma=" << vol << " T=" << maturity << "\n\n";

        std::cout << "== European options ==\n";
        std::cout << std::left << std::setw(8) << "type" << std::setw(14) << "BS price" << std::setw(14)
                  << "MC price" << std::setw(12) << "MC s.e." << "\n";
        for (OptionType type : {OptionType::Call, OptionType::Put}) {
            const VanillaOption o{type, strike, maturity};
            const double exact = bs_price(market, o, vol);
            const McEstimate mc = mc_european(market, o, vol, cfg);
            std::cout << std::setw(8) << (type == OptionType::Call ? "call" : "put") << std::setw(14) << exact
                      << std::setw(14) << mc.price << std::setw(12) << mc.std_error << "\n";
        }

        std::cout << "\n== Greeks (analytic | finite difference) ==\n";
        std::cout << std::setw(8) << "type" << std::setw(24) << "delta" << std::setw(24) << "gamma"
                  << std::setw(24) << "vega" << std::setw(24) << "theta" << std::setw(24) << "rho" << "\n";
        for (OptionType type : {OptionType::Call, OptionType::Put}) {
            const VanillaOption o{type, strike, maturity};
            const Greeks a = bs_greeks(market, o, vol);
            const Greeks f = fd_greeks_vanilla(market, o, vol);
            auto cell = [](double x, double y) {
                std::ostringstream os;
                os << std::setprecision(6) << std::fixed << x << " | " << y;
                return os.str();
            };
            std::cout << std::setw(8) << (type == OptionType::Call ? "call" : "put") << std::setw(24)
                      << cell(a.delta, f.delta) << std::setw(24) << cell(a.gamma, f.gamma) << std::setw(24)
                      << cell(a.vega, f.vega) << std::setw(24) << cell(a.theta, f.theta) << std::setw(24)
                      << cell(a.rho, f.rho) << "\n";
        }

        std::cout << "\n== Barrier options, strike " << strike << ", barrier " << barrier << " ==\n";
        std::cout << std::setw(8) << "type" << std::setw(14) << "closed form" << std::setw(14) << "MC price"
                  << std::setw(12) << "MC s.e." << "\n";
        const BarrierType types[] = {BarrierType::DownAndOut, BarrierType::DownAndIn};
        for (BarrierType bt : types) {
            for (OptionType type : {OptionType::Call, OptionType::Put}) {
                const BarrierOption o{type, bt, strike, barrier, maturity};
                const double exact = barrier_closed_form(market, o, vol);
                const McEstimate mc = mc_barrier(market, o, vol, cfg);
                const std::string name = std::string(bt == BarrierType::DownAndOut ? "DO " : "DI ") +
                                         (type == OptionType::Call ? "call" : "put");
                std::cout << std::setw(8) << name << std::setw(14) << exact << std::setw(14) << mc.price
                          << std::setw(12) << mc.std_error << "\n";
            }
        }
        const double vanilla = bs_price(market, {OptionType::Call, strike, maturity}, vol);
        const double do_call = barrier_closed_form(market, {OptionType::Call, BarrierType::DownAndOut, strike, barrier, maturity}, vol);
        const double di_call = barrier_closed_form(market, {OptionType::Call, BarrierType::DownAndIn, strike, barrier, maturity}, vol);
        std::cout << "\nIn-out parity check (call): DO + DI - vanilla = " << std::scientific << std::setprecision(2)
                  << (do_call + di_call - vanilla) << "\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << "\n";
        print_usage();
        return 1;
    }
}
