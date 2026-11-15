#include "common.hpp"
#include <omp.h>

int main(int argc, char* argv[]) {
    int m = 1024, n = 1024, p = 1024;
    int num_threads = omp_get_max_threads();
    bool json_output = false;

    if (argc >= 4) {
        m = std::atoi(argv[1]);
        n = std::atoi(argv[2]);
        p = std::atoi(argv[3]);
    }
    if (argc >= 5 && argv[4][0] != '-') {
        num_threads = std::atoi(argv[4]);
    }
    for (int i = 1; i < argc; ++i) {
        if (std::string(argv[i]) == "--json") json_output = true;
    }

    std::vector<double> A(static_cast<size_t>(m) * n);
    std::vector<double> B(static_cast<size_t>(n) * p);
    std::vector<double> C(static_cast<size_t>(m) * p, 0.0);

    init_matrix_A(A.data(), m, n);
    init_matrix_B(B.data(), n, p);

    double t0 = omp_get_wtime();
    #pragma omp parallel for schedule(static) num_threads(num_threads)
    for (int i = 0; i < m; ++i) {
        for (int k = 0; k < n; ++k) {
            double a_ik = A[i * n + k];
            for (int j = 0; j < p; ++j) {
                C[i * p + j] += a_ik * B[k * p + j];
            }
        }
    }
    double elapsed_sec = omp_get_wtime() - t0;

    double gflops = calculate_gflops(m, n, p, elapsed_sec);
    double checksum = compute_checksum(C.data(), m, p);

    if (json_output) {
        std::cout << "{\"method\": \"openmp\", \"m\": " << m << ", \"n\": " << n << ", \"p\": " << p 
                  << ", \"threads\": " << num_threads << ", \"time_sec\": " << std::fixed << std::setprecision(6) << elapsed_sec 
                  << ", \"gflops\": " << std::setprecision(4) << gflops 
                  << ", \"checksum\": " << std::setprecision(4) << checksum << "}\n";
    } else {
        std::cout << "OpenMP (" << num_threads << "T): " << elapsed_sec << "s (" << gflops << " GFLOPS)\n";
    }
    return 0;
}
