#ifndef COMMON_HPP
#define COMMON_HPP

#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <chrono>
#include <string>
#include <cstring>
#include <algorithm>

inline void init_matrix_A(double* A, int m, int n) {
    for (int i = 0; i < m; ++i) {
        for (int k = 0; k < n; ++k) {
            A[i * n + k] = static_cast<double>((i + 2 * k) % 100) * 0.01 + 0.05;
        }
    }
}

inline void init_matrix_B(double* B, int n, int p) {
    for (int k = 0; k < n; ++k) {
        for (int j = 0; j < p; ++j) {
            B[k * p + j] = static_cast<double>((2 * k - j + 200) % 100) * 0.01 + 0.05;
        }
    }
}

inline void zero_matrix(double* C, int m, int p) {
    std::memset(C, 0, sizeof(double) * m * p);
}

inline double compute_checksum(const double* C, int m, int p) {
    double sum = 0.0;
    for (int i = 0; i < m * p; ++i) {
        sum += C[i];
    }
    return sum;
}

inline double calculate_gflops(int m, int n, int p, double elapsed_sec) {
    if (elapsed_sec <= 0.0) return 0.0;
    double total_flops = 2.0 * static_cast<double>(m) * static_cast<double>(n) * static_cast<double>(p);
    return (total_flops / elapsed_sec) / 1e9;
}

inline void matmul_ikj_serial(const double* A, const double* B, double* C, int m, int n, int p) {
    zero_matrix(C, m, p);
    for (int i = 0; i < m; ++i) {
        for (int k = 0; k < n; ++k) {
            double a_ik = A[i * n + k];
            for (int j = 0; j < p; ++j) {
                C[i * p + j] += a_ik * B[k * p + j];
            }
        }
    }
}

inline bool verify_matrices(const double* C_ref, const double* C_test, int m, int p, double tol = 1e-6) {
    for (int i = 0; i < m * p; ++i) {
        if (std::abs(C_ref[i] - C_test[i]) > tol) {
            return false;
        }
    }
    return true;
}

class Timer {
public:
    void start() {
        t_start = std::chrono::high_resolution_clock::now();
    }
    double stop() {
        auto t_end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> diff = t_end - t_start;
        return diff.count();
    }
private:
    std::chrono::time_point<std::chrono::high_resolution_clock> t_start;
};

#endif // COMMON_HPP
