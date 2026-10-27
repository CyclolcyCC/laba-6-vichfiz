/**
 * =====================================================================================
 *  BMSTU | Department of Physics (FN-4) | Computational Physics
 *  Laboratory Work No. 6: "Monte Carlo Integration and Variance Reduction"
 *  Module: monte_carlo.hpp
 * =====================================================================================
 */

#ifndef MONTE_CARLO_HPP
#define MONTE_CARLO_HPP

#include <iostream>
#include <iomanip>
#include <cmath>
#include <functional>
#include <string>
#include <vector>
#include <random>
#include <chrono>
#include <sstream>
#include <algorithm>
#include <numbers>

namespace mc {

using MathFunc = std::function<double(double)>;

// -------------------------------------------------------------------------------------
// 1. NUMERICALLY STABLE SINGLE-PASS ACCUMULATOR (WELFORD'S ALGORITHM)
// -------------------------------------------------------------------------------------
struct WelfordAccumulator {
    uint64_t count = 0;
    double mean = 0.0;
    double M2 = 0.0; // Sum of squared differences from current mean

    void update(double x) {
        count++;
        double delta = x - mean;
        mean += delta / count;
        double delta2 = x - mean;
        M2 += delta * delta2;
    }

    [[nodiscard]] double get_mean() const { 
        return mean; 
    }

    [[nodiscard]] double get_variance() const { 
        return (count > 1) ? (M2 / (count - 1)) : 0.0; 
    }

    [[nodiscard]] double get_std_dev() const { 
        return std::sqrt(get_variance()); 
    }

    [[nodiscard]] double get_std_error() const { 
        return (count > 0) ? std::sqrt(get_variance() / count) : 0.0; 
    }

    [[nodiscard]] double get_ci95_halfwidth() const { 
        // 95% two-sided normal quantile: 1.95996398454
        return 1.95996398454 * get_std_error(); 
    }

    [[nodiscard]] double get_ci997_halfwidth() const { 
        // 3-sigma (99.73%) normal quantile: 3.0
        return 3.0 * get_std_error(); 
    }
};

// -------------------------------------------------------------------------------------
// 2. MONTE CARLO INTEGRATION RESULT STRUCT
// -------------------------------------------------------------------------------------
struct MCResult {
    std::string method_name;
    uint64_t n_evals = 0;        // Total function evaluations
    double estimate = 0.0;       // Estimated integral value
    double variance = 0.0;       // Sample variance (per-estimate or normalized)
    double std_error = 0.0;      // Standard error of the mean
    double ci95_lower = 0.0;     // 95% Confidence Interval lower bound
    double ci95_upper = 0.0;     // 95% Confidence Interval upper bound
    double exact_value = 0.0;    // Analytical benchmark
    double abs_error = 0.0;      // |estimate - exact|
    double rel_error = 0.0;      // abs_error / |exact|
    double elapsed_ms = 0.0;     // Execution time in milliseconds
};

// -------------------------------------------------------------------------------------
// 3. FAST 64-BIT RANDOM NUMBER GENERATOR WRAPPER
// -------------------------------------------------------------------------------------
class RandomGenerator {
private:
    std::mt19937_64 rng_;
    std::uniform_real_distribution<double> dist01_{0.0, 1.0};

public:
    explicit RandomGenerator(uint64_t seed = 123456789ULL) : rng_(seed) {}

    void set_seed(uint64_t seed) { 
        rng_.seed(seed); 
    }

    void seed_random_device() {
        std::random_device rd;
        rng_.seed((static_cast<uint64_t>(rd()) << 32) ^ rd());
    }

    [[nodiscard]] double uniform01() { 
        return dist01_(rng_); 
    }

    [[nodiscard]] double uniform(double a, double b) { 
        return a + (b - a) * dist01_(rng_); 
    }
};

// -------------------------------------------------------------------------------------
// 4. MONTE CARLO INTEGRATION METHODS
// -------------------------------------------------------------------------------------

/**
 * @brief 1. Hit-or-Miss (Geometric) Monte Carlo Method
 * Encloses the function under y in [0, M], throws N uniform 2D points into [a,b] x [0,M].
 */
inline MCResult mc_hit_or_miss(const MathFunc& f, double a, double b, double M,
                               uint64_t n_samples, double exact_val, RandomGenerator& rng) {
    auto t_start = std::chrono::high_resolution_clock::now();
    const double area_box = M * (b - a);

    WelfordAccumulator acc;
    for (uint64_t i = 0; i < n_samples; ++i) {
        double x = a + (b - a) * rng.uniform01();
        double y = M * rng.uniform01();
        double hit_val = (y <= f(x)) ? area_box : 0.0;
        acc.update(hit_val);
    }

    auto t_end = std::chrono::high_resolution_clock::now();
    double ms = std::chrono::duration<double, std::milli>(t_end - t_start).count();

    MCResult res;
    res.method_name = "Hit-or-Miss (Geometric)";
    res.n_evals = n_samples;
    res.estimate = acc.get_mean();
    res.variance = acc.get_variance();
    res.std_error = acc.get_std_error();
    res.ci95_lower = res.estimate - acc.get_ci95_halfwidth();
    res.ci95_upper = res.estimate + acc.get_ci95_halfwidth();
    res.exact_value = exact_val;
    res.abs_error = std::abs(res.estimate - exact_val);
    res.rel_error = (exact_val != 0.0) ? (res.abs_error / std::abs(exact_val)) : 0.0;
    res.elapsed_ms = ms;
    return res;
}

/**
 * @brief 2. Sample Mean Monte Carlo Method
 * Computes I = (b - a) * E[f(X)], where X ~ U(a, b).
 */
inline MCResult mc_sample_mean(const MathFunc& f, double a, double b,
                               uint64_t n_samples, double exact_val, RandomGenerator& rng) {
    auto t_start = std::chrono::high_resolution_clock::now();
    const double length = b - a;

    WelfordAccumulator acc;
    for (uint64_t i = 0; i < n_samples; ++i) {
        double x = a + length * rng.uniform01();
        acc.update(length * f(x));
    }

    auto t_end = std::chrono::high_resolution_clock::now();
    double ms = std::chrono::duration<double, std::milli>(t_end - t_start).count();

    MCResult res;
    res.method_name = "Sample Mean (Standard MC)";
    res.n_evals = n_samples;
    res.estimate = acc.get_mean();
    res.variance = acc.get_variance();
    res.std_error = acc.get_std_error();
    res.ci95_lower = res.estimate - acc.get_ci95_halfwidth();
    res.ci95_upper = res.estimate + acc.get_ci95_halfwidth();
    res.exact_value = exact_val;
    res.abs_error = std::abs(res.estimate - exact_val);
    res.rel_error = (exact_val != 0.0) ? (res.abs_error / std::abs(exact_val)) : 0.0;
    res.elapsed_ms = ms;
    return res;
}

/**
 * @brief 3. Antithetic Variates Method
 * Generates negatively correlated pairs: x1 = a + (b-a)*u, x2 = a + (b-a)*(1-u).
 * Uses n_samples/2 pairs (total N function evaluations).
 */
inline MCResult mc_antithetic_variates(const MathFunc& f, double a, double b,
                                       uint64_t n_samples, double exact_val, RandomGenerator& rng) {
    auto t_start = std::chrono::high_resolution_clock::now();
    const double length = b - a;
    const uint64_t n_pairs = n_samples / 2;

    WelfordAccumulator acc;
    for (uint64_t i = 0; i < n_pairs; ++i) {
        double u = rng.uniform01();
        double x1 = a + length * u;
        double x2 = a + length * (1.0 - u);
        double pair_val = 0.5 * length * (f(x1) + f(x2));
        acc.update(pair_val);
    }

    auto t_end = std::chrono::high_resolution_clock::now();
    double ms = std::chrono::duration<double, std::milli>(t_end - t_start).count();

    MCResult res;
    res.method_name = "Antithetic Variates";
    res.n_evals = n_pairs * 2;
    res.estimate = acc.get_mean();
    // Normalized variance per single evaluation: Var(single_eval_equivalent) = 2 * Var(pair)
    res.variance = 2.0 * acc.get_variance();
    res.std_error = acc.get_std_error(); // std_dev(pair) / sqrt(n_pairs)
    res.ci95_lower = res.estimate - acc.get_ci95_halfwidth();
    res.ci95_upper = res.estimate + acc.get_ci95_halfwidth();
    res.exact_value = exact_val;
    res.abs_error = std::abs(res.estimate - exact_val);
    res.rel_error = (exact_val != 0.0) ? (res.abs_error / std::abs(exact_val)) : 0.0;
    res.elapsed_ms = ms;
    return res;
}

/**
 * @brief 4. Control Variates Method
 * Uses known integral I_h = \int_a^b h(x) dx, simulates I_h + (b-a)*[f(X) - h(X)].
 */
inline MCResult mc_control_variates(const MathFunc& f, const MathFunc& h, double I_h,
                                    double a, double b, uint64_t n_samples,
                                    double exact_val, RandomGenerator& rng) {
    auto t_start = std::chrono::high_resolution_clock::now();
    const double length = b - a;

    WelfordAccumulator acc;
    for (uint64_t i = 0; i < n_samples; ++i) {
        double x = a + length * rng.uniform01();
        double diff = f(x) - h(x);
        acc.update(I_h + length * diff);
    }

    auto t_end = std::chrono::high_resolution_clock::now();
    double ms = std::chrono::duration<double, std::milli>(t_end - t_start).count();

    MCResult res;
    res.method_name = "Control Variates";
    res.n_evals = n_samples;
    res.estimate = acc.get_mean();
    res.variance = acc.get_variance();
    res.std_error = acc.get_std_error();
    res.ci95_lower = res.estimate - acc.get_ci95_halfwidth();
    res.ci95_upper = res.estimate + acc.get_ci95_halfwidth();
    res.exact_value = exact_val;
    res.abs_error = std::abs(res.estimate - exact_val);
    res.rel_error = (exact_val != 0.0) ? (res.abs_error / std::abs(exact_val)) : 0.0;
    res.elapsed_ms = ms;
    return res;
}

/**
 * @brief 5. Stratified Sampling Method
 * Subdivides [a, b] into K strata, draws n_samples/K samples in each stratum.
 */
inline MCResult mc_stratified(const MathFunc& f, double a, double b,
                              uint64_t n_samples, uint64_t k_strata,
                              double exact_val, RandomGenerator& rng) {
    auto t_start = std::chrono::high_resolution_clock::now();
    if (k_strata == 0) k_strata = 1;
    uint64_t samples_per_stratum = std::max<uint64_t>(1, n_samples / k_strata);
    uint64_t total_samples = samples_per_stratum * k_strata;

    const double h_stratum = (b - a) / k_strata;
    double total_estimate = 0.0;
    double total_var_of_sum = 0.0;

    for (uint64_t k = 0; k < k_strata; ++k) {
        double sub_a = a + k * h_stratum;
        WelfordAccumulator sub_acc;
        for (uint64_t i = 0; i < samples_per_stratum; ++i) {
            double x = sub_a + h_stratum * rng.uniform01();
            sub_acc.update(h_stratum * f(x));
        }
        total_estimate += sub_acc.get_mean();
        total_var_of_sum += sub_acc.get_variance() / samples_per_stratum;
    }

    auto t_end = std::chrono::high_resolution_clock::now();
    double ms = std::chrono::duration<double, std::milli>(t_end - t_start).count();

    MCResult res;
    res.method_name = "Stratified Sampling";
    res.n_evals = total_samples;
    res.estimate = total_estimate;
    res.std_error = std::sqrt(total_var_of_sum);
    // Effective normalized single-draw variance: Var_eff = total_samples * Var(mean)
    res.variance = total_samples * total_var_of_sum;
    res.ci95_lower = res.estimate - 1.95996398454 * res.std_error;
    res.ci95_upper = res.estimate + 1.95996398454 * res.std_error;
    res.exact_value = exact_val;
    res.abs_error = std::abs(res.estimate - exact_val);
    res.rel_error = (exact_val != 0.0) ? (res.abs_error / std::abs(exact_val)) : 0.0;
    res.elapsed_ms = ms;
    return res;
}

// -------------------------------------------------------------------------------------
// 5. CLASSICAL NEWTON-COTES QUADRATURE RULES (FROM LAB 5 BENCHMARK)
// -------------------------------------------------------------------------------------

inline double integrate_rect_left(const MathFunc& f, double a, double b, int n) {
    const double h = (b - a) / n;
    double sum = 0.0;
    for (int i = 0; i < n; ++i) sum += f(a + i * h);
    return sum * h;
}

inline double integrate_rect_mid(const MathFunc& f, double a, double b, int n) {
    const double h = (b - a) / n;
    double sum = 0.0;
    for (int i = 0; i < n; ++i) sum += f(a + (i + 0.5) * h);
    return sum * h;
}

inline double integrate_trapezoid(const MathFunc& f, double a, double b, int n) {
    const double h = (b - a) / n;
    double sum = 0.5 * (f(a) + f(b));
    for (int i = 1; i < n; ++i) sum += f(a + i * h);
    return sum * h;
}

inline double integrate_simpson(const MathFunc& f, double a, double b, int n) {
    if (n % 2 != 0) n += 1;
    const double h = (b - a) / n;
    double sum = f(a) + f(b);
    for (int i = 1; i < n; ++i) {
        double x = a + i * h;
        sum += (i % 2 == 1) ? 4.0 * f(x) : 2.0 * f(x);
    }
    return (h / 3.0) * sum;
}

} // namespace mc

#endif // MONTE_CARLO_HPP
