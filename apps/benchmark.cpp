// Accuracy and performance benchmark.
//
//   1. European Monte Carlo error versus Black-Scholes for several path counts, with each
//      variance-reduction combination (plain, antithetic, control variate, both).
//   2. Barrier Monte Carlo error versus the closed form for all eight barrier cells.
//   3. Runtime scaling with path count and with thread count.
//
// Writes accuracy_european.csv, accuracy_barrier.csv and scaling.csv to the output directory.
//
// Usage: benchmark [output_dir]

#include "optpricing/barrier.hpp"
#include "optpricing/black_scholes.hpp"
#include "optpricing/monte_carlo.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

using namespace optpricing;

namespace {

using Clock = std::chrono::steady_clock;

MonteCarloConfig make_config(std::size_t paths, bool anti, bool cv, unsigned threads, std::size_t steps = 252) {
    MonteCarloConfig cfg;
    cfg.paths = paths;
    cfg.antithetic = anti;
    cfg.control_variate = cv;
    cfg.threads = threads;
    cfg.steps = steps;
    cfg.seed = 20260101;
    return cfg;
}

void write_csv(const std::filesystem::path& file, const std::string& header, const std::vector<std::string>& rows) {
    std::ofstream out(file);
    if (!out) {
        throw std::runtime_error("cannot open " + file.string() + " for writing");
    }
    out << header << "\n";
    for (const auto& r : rows) {
        out << r << "\n";
    }
}

std::string num(double x, int precision = 8) {
    std::ostringstream os;
    os << std::setprecision(precision) << x;
    return os.str();
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const std::filesystem::path out_dir = argc > 1 ? argv[1] : "benchmark_out";
        std::filesystem::create_directories(out_dir);

        const Market market{100.0, 0.05, 0.02};
        const double vol = 0.25;
        const double T = 1.0;
        const unsigned hw = std::max(1u, std::thread::hardware_concurrency());

        std::cout << "Hardware threads: " << hw << "\n\n";

        // ---- 1. European accuracy --------------------------------------------------------------
        std::cout << "== 1. European Monte Carlo accuracy (S=100, K=100, sigma=0.25, T=1) ==\n";
        std::cout << std::left << std::setw(8) << "type" << std::setw(10) << "paths" << std::setw(12) << "variant"
                  << std::setw(14) << "abs error" << std::setw(14) << "std error" << std::setw(16) << "s.e. ratio"
                  << "\n";
        std::vector<std::string> euro_rows;
        for (OptionType type : {OptionType::Call, OptionType::Put}) {
            const VanillaOption o{type, 100.0, T};
            const double exact = bs_price(market, o, vol);
            const char* tname = type == OptionType::Call ? "call" : "put";
            for (std::size_t paths : {10000u, 100000u, 1000000u}) {
                const double plain_se = mc_european(market, o, vol, make_config(paths, false, false, 1)).std_error;
                struct Variant {
                    const char* name;
                    bool anti;
                    bool cv;
                };
                const Variant variants[] = {{"plain", false, false},
                                            {"antithetic", true, false},
                                            {"control", false, true},
                                            {"both", true, true}};
                for (const Variant& v : variants) {
                    const McEstimate mc = mc_european(market, o, vol, make_config(paths, v.anti, v.cv, 1));
                    const double ratio = plain_se / mc.std_error;
                    std::cout << std::setw(8) << tname << std::setw(10) << paths << std::setw(12) << v.name
                              << std::setw(14) << std::scientific << std::setprecision(3) << std::fabs(mc.price - exact)
                              << std::setw(14) << mc.std_error << std::setw(16) << std::fixed << std::setprecision(2)
                              << ratio << "\n";
                    euro_rows.push_back(std::string(tname) + "," + std::to_string(paths) + "," + v.name + "," +
                                        num(mc.price) + "," + num(exact) + "," + num(std::fabs(mc.price - exact)) +
                                        "," + num(mc.std_error));
                }
            }
        }
        write_csv(out_dir / "accuracy_european.csv", "type,paths,variant,mc_price,bs_price,abs_error,std_error",
                  euro_rows);

        // ---- 2. Barrier accuracy ---------------------------------------------------------------
        std::cout << "\n== 2. Barrier Monte Carlo vs closed form (K=100, H in {85, 115}, 252 steps, 200k paths) ==\n";
        std::cout << std::left << std::setw(8) << "cell" << std::setw(8) << "type" << std::setw(14) << "closed form"
                  << std::setw(14) << "MC price" << std::setw(12) << "std error" << std::setw(10) << "z-score"
                  << "\n";
        struct BCell {
            BarrierType bt;
            double barrier;
        };
        const BCell bcells[] = {{BarrierType::DownAndOut, 85.0}, {BarrierType::DownAndIn, 85.0},
                                {BarrierType::UpAndOut, 115.0},  {BarrierType::UpAndIn, 115.0}};
        std::vector<std::string> barrier_rows;
        for (const BCell& bc : bcells) {
            for (OptionType type : {OptionType::Call, OptionType::Put}) {
                const BarrierOption o{type, bc.bt, 100.0, bc.barrier, T};
                const double exact = barrier_closed_form(market, o, vol);
                const McEstimate mc = mc_barrier(market, o, vol, make_config(200000, true, true, 0));
                const double z = mc.std_error > 0.0 ? (mc.price - exact) / mc.std_error : 0.0;
                const std::string name = std::string(bc.bt == BarrierType::DownAndOut   ? "DO"
                                                     : bc.bt == BarrierType::DownAndIn ? "DI"
                                                     : bc.bt == BarrierType::UpAndOut  ? "UO"
                                                                                       : "UI");
                const char* tname = type == OptionType::Call ? "call" : "put";
                std::cout << std::setw(8) << name << std::setw(8) << tname << std::setw(14) << std::fixed
                          << std::setprecision(6) << exact << std::setw(14) << mc.price << std::setw(12)
                          << std::setprecision(6) << mc.std_error << std::setw(10) << std::setprecision(2) << z << "\n";
                barrier_rows.push_back(name + "," + tname + "," + num(exact) + "," + num(mc.price) + "," +
                                       num(mc.std_error) + "," + num(z));
            }
        }
        write_csv(out_dir / "accuracy_barrier.csv", "cell,type,closed_form,mc_price,std_error,z_score", barrier_rows);

        // ---- 3. Runtime scaling ----------------------------------------------------------------
        std::cout << "\n== 3. Runtime scaling ==\n";
        std::vector<std::string> scale_rows;
        const VanillaOption call{OptionType::Call, 100.0, T};
        const BarrierOption dout{OptionType::Call, BarrierType::DownAndOut, 100.0, 85.0, T};

        std::cout << "European, paths scaling (1 thread):\n";
        for (std::size_t paths : {250000u, 1000000u, 4000000u}) {
            const auto t0 = Clock::now();
            (void)mc_european(market, call, vol, make_config(paths, true, true, 1));
            const double s = std::chrono::duration<double>(Clock::now() - t0).count();
            std::cout << "  paths=" << std::setw(9) << paths << "  time=" << std::fixed << std::setprecision(4) << s
                      << " s  throughput=" << std::setprecision(2) << (static_cast<double>(paths) / s / 1e6)
                      << " Mpaths/s\n";
            scale_rows.push_back("european_paths," + std::to_string(paths) + ",1," + num(s, 6));
        }

        std::cout << "European, thread scaling (4M paths):\n";
        double e1 = 0.0;
        for (unsigned threads : {1u, 2u, 4u, 8u}) {
            const auto t0 = Clock::now();
            (void)mc_european(market, call, vol, make_config(4000000, true, true, threads));
            const double s = std::chrono::duration<double>(Clock::now() - t0).count();
            if (threads == 1) {
                e1 = s;
            }
            std::cout << "  threads=" << threads << "  time=" << std::fixed << std::setprecision(4) << s
                      << " s  speedup=" << std::setprecision(2) << (e1 / s) << "x  efficiency="
                      << std::setprecision(0) << 100.0 * (e1 / s) / threads << "%\n";
            scale_rows.push_back("european_threads,4000000," + std::to_string(threads) + "," + num(s, 6));
        }

        std::cout << "Barrier (252 steps), thread scaling, 400k paths:\n";
        double t1 = 0.0;
        for (unsigned threads : {1u, 2u, 4u, 8u}) {
            const auto t0 = Clock::now();
            (void)mc_barrier(market, dout, vol, make_config(400000, true, true, threads));
            const double s = std::chrono::duration<double>(Clock::now() - t0).count();
            if (threads == 1) {
                t1 = s;
            }
            std::cout << "  threads=" << threads << "  time=" << std::fixed << std::setprecision(4) << s
                      << " s  speedup=" << std::setprecision(2) << (t1 / s) << "x  efficiency="
                      << std::setprecision(0) << 100.0 * (t1 / s) / threads << "%\n";
            scale_rows.push_back("barrier_threads," + std::to_string(threads) + "," + std::to_string(threads) + "," +
                                 num(s, 6));
        }
        write_csv(out_dir / "scaling.csv", "experiment,paths_or_threads,threads,seconds", scale_rows);

        std::cout << "\nCSV results written to " << out_dir.string() << "\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << "\n";
        return 1;
    }
}
