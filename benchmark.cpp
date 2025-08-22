#include <iostream>
#include <vector>
#include <thread>
#include <atomic>
#include <chrono>
#include <random>
#include <algorithm>
#include <iomanip>
#include <cmath>
#include <cstdint>
#include <cstring> // for strcmp
#include <string>

// Parameters
constexpr int64_t K = 10;
constexpr int64_t VECTOR_SIZE = 1'000'00 * K;
constexpr double BASELINE_ITERATIONS = 250.0;
constexpr int64_t X = 42;
constexpr int64_t X_SIZE = 10;

// CRC64-like computation
int64_t compute_crc64(int64_t a, int64_t b, int64_t c, int64_t d) {
    int64_t crc = a ^ b ^ c ^ d;
    for (int64_t i = 0; i < 64; ++i)
        crc = (crc >> 1) ^ (-static_cast<int64_t>(crc & 1) & 0xC96C5795D7870F42LL);
    return crc;
}

// Logarithmic + sqrt math computation
double compute_math_log(int64_t a, int64_t b, int64_t c, int64_t d) {
    double x = static_cast<double>(a + 1);
    double y = static_cast<double>(b + 1);
    double z = static_cast<double>(c + d + 1);
    return std::log2(x) * std::sqrt(y) / std::log1p(z);
}

// Heavy combined computation block
int64_t heavy_block(int64_t a, int64_t b, int64_t c, int64_t d) {
    int64_t crc = compute_crc64(a, b, c, d);
    double math = compute_math_log(a, b, c, d);
    volatile double combined = crc * math;
    return static_cast<int64_t>(combined);
}

// Vector modifications
void vector_modify(std::vector<int64_t>& data, int64_t add_val) {
    for (int64_t i = 0; i < static_cast<int64_t>(data.size()); ++i) {
        data[i] += add_val;
        data[i] ^= (data[i] >> 3);
        data[i] *= 0x9E3779B97F4A7C15LL;
    }
}

// Expand vector by factor X_SIZE
void expand_vector_10x(std::vector<int64_t>& data) {
    int64_t original_size = static_cast<int64_t>(data.size());
    data.resize(original_size * X_SIZE);

    for (int64_t repeat = 1; repeat < X_SIZE; ++repeat) {
        int64_t offset = repeat * original_size;
        for (int64_t i = 0; i < original_size; ++i) {
            int64_t val = data[i] + repeat;
            val ^= i;
            data[offset + i] = val;
        }
    }
}

// Worker function executed by each thread
void worker(std::atomic<bool>& stop_flag, std::atomic<int64_t>& counter, std::atomic<int64_t>& sum) {
    std::mt19937_64 rng(std::random_device{}());
    std::uniform_int_distribution<int64_t> dist(1, VECTOR_SIZE);

    while (!stop_flag.load(std::memory_order_relaxed)) {
        std::vector<int64_t> data;
        data.reserve(VECTOR_SIZE);
        for (int64_t i = 0; i < VECTOR_SIZE; ++i) {
            int64_t a = dist(rng);
            int64_t b = dist(rng);
            int64_t c = a + 4;
            int64_t d = b + 11;

            int64_t value = heavy_block(a, b, c, d);
            data.push_back(value);
        }
        std::sort(data.begin(), data.end());
        expand_vector_10x(data);
        vector_modify(data, X);
        counter.fetch_add(1, std::memory_order_relaxed);
        sum.fetch_add(data.at((dist(rng) % VECTOR_SIZE) * X_SIZE), std::memory_order_relaxed);
    }
}

struct BenchmarkResult {
    int64_t total_iterations;
    std::vector<int64_t> iterations_per_thread;
};

// Run single benchmark with given threads and duration
BenchmarkResult run_single_benchmark(int64_t num_threads, int64_t duration_seconds) {
    std::atomic<bool> stop_flag{false};
    std::vector<std::thread> threads;
    std::vector<std::atomic<int64_t>> iterations(num_threads);
    std::vector<std::atomic<int64_t>> sums(num_threads);

    for (auto& c : iterations) c.store(0);
    for (auto& s : sums) s.store(0);

    for (int64_t i = 0; i < num_threads; ++i)
        threads.emplace_back(worker, std::ref(stop_flag), std::ref(iterations[i]), std::ref(sums[i]));

    std::this_thread::sleep_for(std::chrono::seconds(duration_seconds));
    stop_flag.store(true);

    for (auto& t : threads) t.join();

    int64_t total_iterations = 0;
    std::vector<int64_t> iterations_vec;
    iterations_vec.reserve(num_threads);

    for (const auto& it : iterations) {
        int64_t val = it.load();
        total_iterations += val;
        iterations_vec.push_back(val);
    }

    return {total_iterations, iterations_vec};
}

// Print results of single benchmark run
void print_single_run_results(int64_t total_iterations,
                              std::vector<int64_t>& iterations_vec,
                              int64_t duration_seconds,
                              int64_t num_threads,
                              bool show_threads = false) {
    std::sort(iterations_vec.begin(), iterations_vec.end(), std::greater<int64_t>());

    std::cout << "Total iterations across 1 thread: "
              << total_iterations / num_threads << "\n";
    std::cout << "Total iterations across all threads: "
              << total_iterations << "\n";

    double performance_ratio = static_cast<double>(total_iterations)
                               * (10.0 / duration_seconds)
                               / BASELINE_ITERATIONS;

    std::cout << std::fixed << std::setprecision(2)
              << "Performance: x" << performance_ratio
              << " of 4 vCPU AMD Genoa(zen4)\n";

    if (show_threads) {
        std::cout << "Iterations per thread: ";
        for (auto v : iterations_vec) std::cout << v << " ";
        std::cout << "\n";
    }
}

// Print drop percentages compared to previous
void print_drop_percentages(const std::vector<int64_t>& values) {
    if (values.empty()) return;
    std::cout << "Drop per run vs previous: ";
    int64_t prev = values.front();
    std::cout << "0% ";
    for (size_t i = 1; i < values.size(); ++i) {
        int64_t drop = static_cast<int64_t>((static_cast<double>(prev - values[i]) / prev) * 100 + 0.5);
        std::cout << drop << "% ";
        prev = values[i];
    }
    std::cout << "\n";
}

// Run thread-scaling benchmark
void run_thread_scaling_benchmark(int64_t max_threads, int64_t duration_seconds) {
    std::vector<int64_t> averages_int;
    std::vector<int64_t> total_iters_list;

    for (int64_t t = max_threads; t >= 1; --t) {
        auto result = run_single_benchmark(t, duration_seconds);
        int64_t avg_int = static_cast<int64_t>(result.total_iterations / t + 0.5);
        averages_int.push_back(avg_int);
        total_iters_list.push_back(result.total_iterations);
    }

    std::cout << "Iterations per thread: ";
    for (auto v : averages_int) std::cout << v << " ";
    std::cout << "\n";

    std::cout << "Iterations per run: ";
    for (auto v : total_iters_list) std::cout << v << " ";
    std::cout << "\n";

    print_drop_percentages(total_iters_list);
}

int main(int argc, char* argv[]) {
    int64_t duration_seconds = 10;
    bool is_duration_default = true;
    int64_t num_threads = 0;
    bool thread_scaling_mode = false;
    bool show_threads = false;

    std::cout << "Run with --help to see available options.\n\n";

    for (int i = 1; i < argc; ++i) {
        if (std::strcmp(argv[i], "-s") == 0 && i + 1 < argc)
            duration_seconds = std::stoll(argv[++i]), is_duration_default = false;
        else if (std::strcmp(argv[i], "-th") == 0 && i + 1 < argc)
            num_threads = std::stoll(argv[++i]);
        else if (std::strcmp(argv[i], "--thread-scaling") == 0 || std::strcmp(argv[i], "-ts") == 0)
            thread_scaling_mode = true;
        else if (std::strcmp(argv[i], "--show-threads") == 0 || std::strcmp(argv[i], "-st") == 0)
            show_threads = true;
        else if (std::strcmp(argv[i], "--help") == 0) {
            std::cout << "Usage: ./benchmark [-s seconds] [-th threads] [--thread-scaling] [-ts] [--show-threads] [-st]\n";
            return 0;
        }
    }

    int64_t hw_threads = static_cast<int64_t>(std::thread::hardware_concurrency());
    if (num_threads == 0) num_threads = hw_threads;

    if (!thread_scaling_mode) {
        std::cout << "Benchmark duration: " << duration_seconds << " seconds\n";
        std::cout << "Threads used: " << num_threads << " (hardware: " << hw_threads << ")\n";
        auto result = run_single_benchmark(num_threads, duration_seconds);
        print_single_run_results(result.total_iterations,
                                 result.iterations_per_thread,
                                 duration_seconds,
                                 num_threads,
                                 show_threads);
    } else {
        if (is_duration_default) duration_seconds = 10;
        if (num_threads == 0) num_threads = hw_threads;
        std::cout << "Running thread scaling benchmark for " << duration_seconds
                  << " seconds per run, from " << num_threads << " threads down to 1.\n";
        run_thread_scaling_benchmark(num_threads, duration_seconds);
    }

    return 0;
}
