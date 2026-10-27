/**
 * =====================================================================================
 *  BMSTU | Department of Physics (FN-4) | Computational Physics
 *  Laboratory Work No. 6: "Monte Carlo Integration and Variance Reduction"
 *  Variant 7 Benchmark Comparison (Integral from Lab 5)
 * =====================================================================================
 */

#include "monte_carlo.hpp"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>

namespace {

void print_separator(int len = 125, char ch = '=') {
    std::cout << std::string(len, ch) << "\n";
}

std::string format_sci(double val, int prec = 3) {
    std::ostringstream ss;
    ss << std::scientific << std::setprecision(prec) << val;
    return ss.str();
}

std::string format_fix(double val, int prec = 8) {
    std::ostringstream ss;
    ss << std::fixed << std::setprecision(prec) << val;
    return ss.str();
}

} // anonymous namespace

int main() {
    // Synchronize C++ streams with fast I/O
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(nullptr);

    mc::RandomGenerator rng(42ULL); // Fixed seed for reproducible benchmarks

    // 1. Unified concise CSV for instructor demonstration
    std::ofstream csv_main("lab6_results.csv");
    csv_main << "Категория,Интеграл,Метод,Число_вызовов_N,Оценка_интеграла,Точное_значение,Абс_погрешность,Дисперсия_s2,Снижение_дисперсии\n";

    // 2. Focused table 1: Task 1 (Common integral exp(x))
    std::ofstream csv_t1("table1_monte_carlo.csv");
    csv_t1 << "Метод,Число_вызовов_N,Оценка_интеграла,Точное_значение,Абс_погрешность,Дисперсия_s2,Снижение_дисперсии\n";

    // 3. Focused table 2: Task 2 (Lab 5 Variant 7 comparison)
    std::ofstream csv_t2("table2_comparison_lab5.csv");
    csv_t2 << "Семейство,Метод,Число_вызовов_N,Оценка_интеграла,Точное_значение,Абс_погрешность,Порядок_точности\n";

    // 4. Focused table 3: Convergence study by N
    std::ofstream csv_t3("table3_convergence.csv");
    csv_t3 << "N,SampleMean_Err,Antithetic_Err,ControlVar_Err,Stratified_Err\n";

    std::cout << "\n";
    print_separator(125, '=');
    std::cout << "  BMSTU | DEPARTMENT OF PHYSICS (FN-4) | COMPUTATIONAL PHYSICS\n";
    std::cout << "  LABORATORY WORK NO. 6: MONTE CARLO INTEGRATION & VARIANCE REDUCTION\n";
    std::cout << "  Comparison with Laboratory Work No. 5 (Variant 7)\n";
    print_separator(125, '=');

    // =================================================================================
    // PART 1: COMMON INTEGRAL 1: \int_0^1 \exp(x) dx
    // =================================================================================
    {
        const double a = 0.0;
        const double b = 1.0;
        const double exact_exp = std::numbers::e - 1.0;
        const double M_exp = std::numbers::e; // max exp(x) on [0, 1]

        auto f_exp = [](double x) { return std::exp(x); };
        // Control variate: linear chord h(x) = 1 + (e - 1)*x
        auto h_exp = [exact_exp](double x) { return 1.0 + exact_exp * x; };
        const double I_h_exp = (std::numbers::e + 1.0) / 2.0;

        std::cout << "\n>>> TASK 1.1: Common Integral 1\n";
        std::cout << "    I = \\int_0^1 \\exp(x) dx = e - 1 = " 
                  << std::fixed << std::setprecision(12) << exact_exp << "\n";
        print_separator(125, '-');

        const uint64_t N_benchmark = 1'000'000;
        std::cout << "Benchmark at N = " << N_benchmark << " function evaluations:\n";
        std::cout << std::left
                  << std::setw(28) << "Method"
                  << std::setw(12) << "N_evals"
                  << std::setw(16) << "Estimate"
                  << std::setw(15) << "Variance (s^2)"
                  << std::setw(14) << "Std. Error"
                  << std::setw(14) << "Abs. Error"
                  << std::setw(14) << "Var. Reduct."
                  << std::setw(12) << "Time (ms)"
                  << "\n";
        print_separator(125, '-');

        std::vector<mc::MCResult> results_exp;
        results_exp.push_back(mc::mc_hit_or_miss(f_exp, a, b, M_exp, N_benchmark, exact_exp, rng));
        results_exp.push_back(mc::mc_sample_mean(f_exp, a, b, N_benchmark, exact_exp, rng));
        results_exp.push_back(mc::mc_antithetic_variates(f_exp, a, b, N_benchmark, exact_exp, rng));
        results_exp.push_back(mc::mc_control_variates(f_exp, h_exp, I_h_exp, a, b, N_benchmark, exact_exp, rng));
        results_exp.push_back(mc::mc_stratified(f_exp, a, b, N_benchmark, 100, exact_exp, rng));

        const double base_var_exp = results_exp[1].variance; // Sample mean variance

        for (const auto& res : results_exp) {
            double var_ratio = (res.variance > 0.0) ? (base_var_exp / res.variance) : 1.0;
            std::string ratio_str = (var_ratio < 0.99) ? (format_fix(var_ratio, 2) + "x (хуже)")
                                  : (var_ratio > 1.01) ? (format_fix(var_ratio, 2) + "x")
                                  : "1.00x (база)";

            std::cout << std::left
                      << std::setw(28) << res.method_name
                      << std::setw(12) << res.n_evals
                      << std::setw(16) << std::fixed << std::setprecision(9) << res.estimate
                      << std::setw(15) << format_sci(res.variance, 4)
                      << std::setw(14) << format_sci(res.std_error, 3)
                      << std::setw(14) << format_sci(res.abs_error, 3)
                      << std::setw(14) << std::fixed << std::setprecision(2) << var_ratio << "x"
                      << std::setw(12) << std::fixed << std::setprecision(2) << res.elapsed_ms
                      << "\n";

            // Main CSV
            csv_main << "Задание 1,int_0^1 exp(x) dx,\"" << res.method_name << "\","
                     << res.n_evals << ","
                     << format_fix(res.estimate, 8) << ","
                     << format_fix(exact_exp, 8) << ","
                     << format_sci(res.abs_error, 3) << ","
                     << format_sci(res.variance, 3) << ",\""
                     << ratio_str << "\"\n";

            // Table 1 CSV
            csv_t1 << "\"" << res.method_name << "\","
                   << res.n_evals << ","
                   << format_fix(res.estimate, 8) << ","
                   << format_fix(exact_exp, 8) << ","
                   << format_sci(res.abs_error, 3) << ","
                   << format_sci(res.variance, 3) << ",\""
                   << ratio_str << "\"\n";
        }

        // Convergence study for exp(x)
        std::cout << "\nConvergence study across sample sizes N for \\int_0^1 \\exp(x) dx:\n";
        std::vector<uint64_t> sample_sizes = {100, 1'000, 10'000, 100'000, 1'000'000, 10'000'000};
        
        std::cout << std::left
                  << std::setw(12) << "N"
                  << std::setw(22) << "SampleMean Err"
                  << std::setw(22) << "Antithetic Err"
                  << std::setw(22) << "ControlVar Err"
                  << std::setw(22) << "Stratified Err"
                  << "\n";
        print_separator(90, '-');

        for (uint64_t n : sample_sizes) {
            auto r_sm   = mc::mc_sample_mean(f_exp, a, b, n, exact_exp, rng);
            auto r_anti = mc::mc_antithetic_variates(f_exp, a, b, n, exact_exp, rng);
            auto r_cv   = mc::mc_control_variates(f_exp, h_exp, I_h_exp, a, b, n, exact_exp, rng);
            auto r_str  = mc::mc_stratified(f_exp, a, b, n, std::min<uint64_t>(100, n / 10), exact_exp, rng);

            std::cout << std::left
                      << std::setw(12) << n
                      << std::setw(22) << format_sci(r_sm.abs_error, 4)
                      << std::setw(22) << format_sci(r_anti.abs_error, 4)
                      << std::setw(22) << format_sci(r_cv.abs_error, 4)
                      << std::setw(22) << format_sci(r_str.abs_error, 4)
                      << "\n";

            // Convergence Table 3 CSV
            csv_t3 << n << ","
                   << format_sci(r_sm.abs_error, 3) << ","
                   << format_sci(r_anti.abs_error, 3) << ","
                   << format_sci(r_cv.abs_error, 3) << ","
                   << format_sci(r_str.abs_error, 3) << "\n";
        }
    }

    // =================================================================================
    // PART 2: COMMON INTEGRAL 2: \int_0^{\pi/2} \sin(x) dx
    // =================================================================================
    {
        const double a = 0.0;
        const double b = std::numbers::pi / 2.0;
        const double exact_sin = 1.0;
        const double M_sin = 1.0;

        auto f_sin = [](double x) { return std::sin(x); };
        // Quadratic control variate: h(x) = (4/pi)*x - (4/pi^2)*x^2, I_h = pi/3
        const double pi = std::numbers::pi;
        auto h_sin = [pi](double x) { return (4.0 / pi) * x - (4.0 / (pi * pi)) * (x * x); };
        const double I_h_sin = pi / 3.0;

        std::cout << "\n>>> TASK 1.2: Common Integral 2\n";
        std::cout << "    I = \\int_0^{\\pi/2} \\sin(x) dx = 1.000000000000\n";
        print_separator(125, '-');

        const uint64_t N_benchmark = 1'000'000;
        std::vector<mc::MCResult> results_sin;
        results_sin.push_back(mc::mc_hit_or_miss(f_sin, a, b, M_sin, N_benchmark, exact_sin, rng));
        results_sin.push_back(mc::mc_sample_mean(f_sin, a, b, N_benchmark, exact_sin, rng));
        results_sin.push_back(mc::mc_antithetic_variates(f_sin, a, b, N_benchmark, exact_sin, rng));
        results_sin.push_back(mc::mc_control_variates(f_sin, h_sin, I_h_sin, a, b, N_benchmark, exact_sin, rng));
        results_sin.push_back(mc::mc_stratified(f_sin, a, b, N_benchmark, 100, exact_sin, rng));

        const double base_var_sin = results_sin[1].variance;

        std::cout << std::left
                  << std::setw(28) << "Method"
                  << std::setw(12) << "N_evals"
                  << std::setw(16) << "Estimate"
                  << std::setw(15) << "Variance (s^2)"
                  << std::setw(14) << "Std. Error"
                  << std::setw(14) << "Abs. Error"
                  << std::setw(14) << "Var. Reduct."
                  << std::setw(12) << "Time (ms)"
                  << "\n";
        print_separator(125, '-');

        for (const auto& res : results_sin) {
            double var_ratio = (res.variance > 0.0) ? (base_var_sin / res.variance) : 1.0;
            std::string ratio_str = (var_ratio < 0.99) ? (format_fix(var_ratio, 2) + "x (хуже)")
                                  : (var_ratio > 1.01) ? (format_fix(var_ratio, 2) + "x")
                                  : "1.00x (база)";

            std::cout << std::left
                      << std::setw(28) << res.method_name
                      << std::setw(12) << res.n_evals
                      << std::setw(16) << std::fixed << std::setprecision(9) << res.estimate
                      << std::setw(15) << format_sci(res.variance, 4)
                      << std::setw(14) << format_sci(res.std_error, 3)
                      << std::setw(14) << format_sci(res.abs_error, 3)
                      << std::setw(14) << std::fixed << std::setprecision(2) << var_ratio << "x"
                      << std::setw(12) << std::fixed << std::setprecision(2) << res.elapsed_ms
                      << "\n";

            csv_main << "Задание 1,int_0^{pi/2} sin(x) dx,\"" << res.method_name << "\","
                     << res.n_evals << ","
                     << format_fix(res.estimate, 8) << ","
                     << format_fix(exact_sin, 8) << ","
                     << format_sci(res.abs_error, 3) << ","
                     << format_sci(res.variance, 3) << ",\""
                     << ratio_str << "\"\n";
        }
    }

    // =================================================================================
    // PART 3: TASK 2: LAB 5 VARIANT 7 INTEGRAL BENCHMARK
    // =================================================================================
    {
        // Selected: I_1 = \int_0^{\sqrt{3}} \sqrt{3 - x^2} dx = 3*pi/4 \approx 2.356194490192345
        const double a = 0.0;
        const double b = std::sqrt(3.0);
        const double exact_circ = 3.0 * std::numbers::pi / 4.0;
        const double M_circ = std::sqrt(3.0);

        auto f_circ = [](double x) {
            double val = 3.0 - x * x;
            return (val > 0.0) ? std::sqrt(val) : 0.0;
        };

        // Linear chord h(x) = sqrt(3) * (1 - x / sqrt(3)), I_h = 1.5
        auto h_circ = [b](double x) { return b * (1.0 - x / b); };
        const double I_h_circ = 1.5;

        std::cout << "\n>>> TASK 2: Lab 5 (Variant 7) Selection Integral\n";
        std::cout << "    I = \\int_0^{\\sqrt{3}} \\sqrt{3 - x^2} dx = 3\\pi/4 = "
                  << std::fixed << std::setprecision(12) << exact_circ << "\n";
        std::cout << "    [Geometric meaning: Quarter of a circle of radius R = \\sqrt{3}]\n";
        std::cout << "    [Note: Derivative f'(x) -> -\\infty at x = \\sqrt{3}, which broke Simpson in Lab 5]\n";
        print_separator(125, '-');

        const uint64_t N_benchmark = 1'000'000;
        std::vector<mc::MCResult> results_circ;
        results_circ.push_back(mc::mc_hit_or_miss(f_circ, a, b, M_circ, N_benchmark, exact_circ, rng));
        results_circ.push_back(mc::mc_sample_mean(f_circ, a, b, N_benchmark, exact_circ, rng));
        results_circ.push_back(mc::mc_antithetic_variates(f_circ, a, b, N_benchmark, exact_circ, rng));
        results_circ.push_back(mc::mc_control_variates(f_circ, h_circ, I_h_circ, a, b, N_benchmark, exact_circ, rng));
        results_circ.push_back(mc::mc_stratified(f_circ, a, b, N_benchmark, 100, exact_circ, rng));

        const double base_var_circ = results_circ[1].variance;

        std::cout << "Monte Carlo Integration Results (N = " << N_benchmark << "):\n";
        std::cout << std::left
                  << std::setw(28) << "Method"
                  << std::setw(12) << "N_evals"
                  << std::setw(16) << "Estimate"
                  << std::setw(15) << "Variance (s^2)"
                  << std::setw(14) << "Std. Error"
                  << std::setw(14) << "Abs. Error"
                  << std::setw(14) << "Var. Reduct."
                  << std::setw(12) << "Time (ms)"
                  << "\n";
        print_separator(125, '-');

        for (const auto& res : results_circ) {
            double var_ratio = (res.variance > 0.0) ? (base_var_circ / res.variance) : 1.0;
            std::string ratio_str = (var_ratio < 0.99) ? (format_fix(var_ratio, 2) + "x (хуже)")
                                  : (var_ratio > 1.01) ? (format_fix(var_ratio, 2) + "x")
                                  : "1.00x (база)";

            std::cout << std::left
                      << std::setw(28) << res.method_name
                      << std::setw(12) << res.n_evals
                      << std::setw(16) << std::fixed << std::setprecision(9) << res.estimate
                      << std::setw(15) << format_sci(res.variance, 4)
                      << std::setw(14) << format_sci(res.std_error, 3)
                      << std::setw(14) << format_sci(res.abs_error, 3)
                      << std::setw(14) << std::fixed << std::setprecision(2) << var_ratio << "x"
                      << std::setw(12) << std::fixed << std::setprecision(2) << res.elapsed_ms
                      << "\n";

            csv_main << "Задание 2 (ММК),int_0^{sqrt3} sqrt(3-x^2) dx,\"" << res.method_name << "\","
                     << res.n_evals << ","
                     << format_fix(res.estimate, 8) << ","
                     << format_fix(exact_circ, 8) << ","
                     << format_sci(res.abs_error, 3) << ","
                     << format_sci(res.variance, 3) << ",\""
                     << ratio_str << "\"\n";
        }

        // -----------------------------------------------------------------------------
        // DIRECT COMPARISON WITH LAB 5 QUADRATURES (RECTANGLES, TRAPEZOIDS, SIMPSON)
        // -----------------------------------------------------------------------------
        std::cout << "\nDIRECT COMPARISON: Monte Carlo vs Lab 5 Newton-Cotes Quadratures\n";
        print_separator(125, '-');
        std::cout << std::left
                  << std::setw(32) << "Integration Family & Method"
                  << std::setw(14) << "Evaluations"
                  << std::setw(18) << "Estimate"
                  << std::setw(16) << "Abs. Error"
                  << std::setw(16) << "Effective Order"
                  << std::setw(16) << "CPU Time (ms)"
                  << "\n";
        print_separator(125, '-');

        std::vector<int> quad_intervals = {10, 100, 1'000, 10'000};
        
        for (int n_int : quad_intervals) {
            auto t0 = std::chrono::high_resolution_clock::now();
            double I_rect = mc::integrate_rect_mid(f_circ, a, b, n_int);
            auto t1 = std::chrono::high_resolution_clock::now();
            double ms_rect = std::chrono::duration<double, std::milli>(t1 - t0).count();

            t0 = std::chrono::high_resolution_clock::now();
            double I_trap = mc::integrate_trapezoid(f_circ, a, b, n_int);
            t1 = std::chrono::high_resolution_clock::now();
            double ms_trap = std::chrono::duration<double, std::milli>(t1 - t0).count();

            t0 = std::chrono::high_resolution_clock::now();
            double I_simp = mc::integrate_simpson(f_circ, a, b, n_int);
            t1 = std::chrono::high_resolution_clock::now();
            double ms_simp = std::chrono::duration<double, std::milli>(t1 - t0).count();

            std::cout << std::left
                      << std::setw(32) << ("Midpoint Rect (N=" + std::to_string(n_int) + ")")
                      << std::setw(14) << n_int
                      << std::setw(18) << std::fixed << std::setprecision(9) << I_rect
                      << std::setw(16) << format_sci(std::abs(I_rect - exact_circ), 4)
                      << std::setw(16) << "O(h^1.5)*"
                      << std::setw(16) << std::fixed << std::setprecision(3) << ms_rect
                      << "\n";

            std::cout << std::left
                      << std::setw(32) << ("Trapezoid (N=" + std::to_string(n_int) + ")")
                      << std::setw(14) << (n_int + 1)
                      << std::setw(18) << std::fixed << std::setprecision(9) << I_trap
                      << std::setw(16) << format_sci(std::abs(I_trap - exact_circ), 4)
                      << std::setw(16) << "O(h^1.5)*"
                      << std::setw(16) << std::fixed << std::setprecision(3) << ms_trap
                      << "\n";

            std::cout << std::left
                      << std::setw(32) << ("Simpson (N=" + std::to_string(n_int) + ")")
                      << std::setw(14) << (n_int + 1)
                      << std::setw(18) << std::fixed << std::setprecision(9) << I_simp
                      << std::setw(16) << format_sci(std::abs(I_simp - exact_circ), 4)
                      << std::setw(16) << "O(h^1.5)* (degraded)"
                      << std::setw(16) << std::fixed << std::setprecision(3) << ms_simp
                      << "\n";

            // Only export N=100 and N=10000 into Table 2 for a crisp, neat comparison
            if (n_int == 100 || n_int == 10000) {
                csv_t2 << "Ньютон-Котес (Лаб 5),\"Средние прямоугольники (N=" << n_int << ")\","
                       << n_int << ","
                       << format_fix(I_rect, 8) << ","
                       << format_fix(exact_circ, 8) << ","
                       << format_sci(std::abs(I_rect - exact_circ), 3) << ",\"O(h^1.5)\"\n";

                csv_t2 << "Ньютон-Котес (Лаб 5),\"Формула трапеций (N=" << n_int << ")\","
                       << (n_int + 1) << ","
                       << format_fix(I_trap, 8) << ","
                       << format_fix(exact_circ, 8) << ","
                       << format_sci(std::abs(I_trap - exact_circ), 3) << ",\"O(h^1.5)\"\n";

                csv_t2 << "Ньютон-Котес (Лаб 5),\"Формула Симпсона (N=" << n_int << ")\","
                       << (n_int + 1) << ","
                       << format_fix(I_simp, 8) << ","
                       << format_fix(exact_circ, 8) << ","
                       << format_sci(std::abs(I_simp - exact_circ), 3) << ",\"O(h^1.5) деградация\"\n";
            }
        }

        print_separator(125, '-');
        auto mc_high = mc::mc_antithetic_variates(f_circ, a, b, 10'000'000, exact_circ, rng);
        auto mc_sm1m = mc::mc_sample_mean(f_circ, a, b, 1'000'000, exact_circ, rng);

        std::cout << std::left
                  << std::setw(32) << "MC Antithetic (N=10^7)"
                  << std::setw(14) << mc_high.n_evals
                  << std::setw(18) << std::fixed << std::setprecision(9) << mc_high.estimate
                  << std::setw(16) << format_sci(mc_high.abs_error, 4)
                  << std::setw(16) << "O(N^-0.5) (exact)"
                  << std::setw(16) << std::fixed << std::setprecision(3) << mc_high.elapsed_ms
                  << "\n";

        // Add Monte Carlo rows into Table 2
        csv_t2 << "Монте-Карло (Лаб 6),\"Выборочное среднее (N=10^6)\","
               << mc_sm1m.n_evals << ","
               << format_fix(mc_sm1m.estimate, 8) << ","
               << format_fix(exact_circ, 8) << ","
               << format_sci(mc_sm1m.abs_error, 3) << ",\"O(N^-0.5) строго\"\n";

        csv_t2 << "Монте-Карло (Лаб 6),\"Антитетические переменные (N=10^7)\","
               << mc_high.n_evals << ","
               << format_fix(mc_high.estimate, 8) << ","
               << format_fix(exact_circ, 8) << ","
               << format_sci(mc_high.abs_error, 3) << ",\"O(N^-0.5) строго\"\n";

        csv_main << "Задание 2 (Квадратуры Лаб 5),int_0^{sqrt3} sqrt(3-x^2) dx,\"Симпсон (N=100)\","
                 << 101 << ","
                 << format_fix(mc::integrate_simpson(f_circ, a, b, 100), 8) << ","
                 << format_fix(exact_circ, 8) << ","
                 << format_sci(std::abs(mc::integrate_simpson(f_circ, a, b, 100) - exact_circ), 3)
                 << ",--,\"Деградация до O(h^1.5)\"\n";

        csv_main << "Задание 2 (Квадратуры Лаб 5),int_0^{sqrt3} sqrt(3-x^2) dx,\"Симпсон (N=10000)\","
                 << 10001 << ","
                 << format_fix(mc::integrate_simpson(f_circ, a, b, 10000), 8) << ","
                 << format_fix(exact_circ, 8) << ","
                 << format_sci(std::abs(mc::integrate_simpson(f_circ, a, b, 10000) - exact_circ), 3)
                 << ",--,\"Деградация до O(h^1.5)\"\n";

        print_separator(125, '=');
        std::cout << " [*] Analysis Note:\n"
                  << "     For f(x) = \\sqrt{3 - x^2}, Simpson's theoretical order O(h^4) dropped to O(h^1.5)\n"
                  << "     due to vertical tangent derivative singularity f'(x) -> -\\infty at x = \\sqrt{3}.\n"
                  << "     Meanwhile, Monte Carlo requires NO derivative continuity and strictly converges\n"
                  << "     at rate O(N^-0.5). However, for 1D integration, Simpson still reaches ~10^-6 error\n"
                  << "     with fewer function calls than basic Monte Carlo.\n";
        print_separator(125, '=');
    }

    csv_main.close();
    csv_t1.close();
    csv_t2.close();
    csv_t3.close();

    std::cout << "\n[SUCCESS] Calculations completed.\n"
              << "  -> 'lab6_results.csv'          : Основная понятная таблица для преподавателя\n"
              << "  -> 'table1_monte_carlo.csv'    : Таблица 1 (ММК методы для exp(x))\n"
              << "  -> 'table2_comparison_lab5.csv': Таблица 2 (Прямое сравнение ММК и Лаб №5 Симпсон/трапеции)\n"
              << "  -> 'table3_convergence.csv'    : Таблица 3 (Сходимость погрешности по N)\n\n";

    return 0;
}
