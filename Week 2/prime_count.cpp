#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <string>
#include <omp.h>

// Kiểm tra số nguyên tố
inline bool is_prime(long long num) {
    if (num <= 1) return false;
    if (num == 2 || num == 3) return true;
    if (num % 2 == 0 || num % 3 == 0) return false;
    for (long long i = 5; i * i <= num; i += 6) {
        if (num % i == 0 || num % (i + 2) == 0) return false;
    }
    return true;
}

// Đếm số nguyên tố song song bằng OpenMP
long long count_primes_parallel(long long n, int num_threads) {
    long long count = 0;
    #pragma omp parallel for reduction(+:count) schedule(guided, 1000) num_threads(num_threads)
    for (long long i = 1; i <= n; ++i) {
        if (is_prime(i)) {
            count++;
        }
    }
    return count;
}

int main(int argc, char* argv[]) {
    long long n = 10000000; // Mặc định 10 triệu
    int num_threads = 4;
    bool json_output = false;

    if (argc >= 2) n = std::atoll(argv[1]);
    if (argc >= 3) num_threads = std::atoi(argv[2]);
    for (int i = 1; i < argc; ++i) {
        if (std::string(argv[i]) == "--json") json_output = true;
    }

    double t0 = omp_get_wtime();
    long long count = count_primes_parallel(n, num_threads);
    double t1 = omp_get_wtime();
    double elapsed_sec = t1 - t0;

    if (json_output) {
        std::cout << "{\"n\": " << n << ", \"threads\": " << num_threads 
                  << ", \"count\": " << count << ", \"time_sec\": " 
                  << std::fixed << std::setprecision(6) << elapsed_sec << "}\n";
    } else {
        std::cout << "========================================================\n";
        std::cout << "  DEM SO NGUYEN TO DEN N (PRIME COUNTING - OPENMP)\n";
        std::cout << "  N = " << n << " | So luong luong OpenMP: " << num_threads << "\n";
        std::cout << "========================================================\n";
        std::cout << "So luong so nguyen to pi(N): " << count << "\n";
        std::cout << "Thoi gian thuc thi: " << std::fixed << std::setprecision(6) 
                  << elapsed_sec << " giay (" << elapsed_sec * 1000.0 << " ms)\n";
        std::cout << "--------------------------------------------------------\n";
    }

    return 0;
}
