#include "common.hpp"
#include <thread>
#include <vector>

void matmul_thread_worker(const double* A, const double* B, double* C, 
                          int r_start, int r_end, int n, int p) {
    for (int i = r_start; i < r_end; ++i) {
        for (int k = 0; k < n; ++k) {
            double a_ik = A[i * n + k];
            for (int j = 0; j < p; ++j) {
                C[i * p + j] += a_ik * B[k * p + j];
            }
        }
    }
}

int main(int argc, char* argv[]) {
    int m = 1024, n = 1024, p = 1024;
    int num_threads = std::thread::hardware_concurrency();
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

    Timer timer;
    timer.start();

    std::vector<std::thread> workers;
    workers.reserve(num_threads);
    int rows_per_thread = m / num_threads;
    int rem = m % num_threads;
    int current_row = 0;

    for (int t = 0; t < num_threads; ++t) {
        int count = rows_per_thread + (t < rem ? 1 : 0);
        int r_start = current_row;
        int r_end = current_row + count;
        current_row = r_end;
        if (count > 0) {
            workers.emplace_back(matmul_thread_worker, A.data(), B.data(), C.data(), r_start, r_end, n, p);
        }
    }
    for (auto& th : workers) {
        if (th.joinable()) th.join();
    }
    double elapsed_sec = timer.stop();

    double gflops = calculate_gflops(m, n, p, elapsed_sec);
    double checksum = compute_checksum(C.data(), m, p);

    if (json_output) {
        std::cout << "{\"method\": \"cpp_thread\", \"m\": " << m << ", \"n\": " << n << ", \"p\": " << p 
                  << ", \"threads\": " << num_threads << ", \"time_sec\": " << std::fixed << std::setprecision(6) << elapsed_sec 
                  << ", \"gflops\": " << std::setprecision(4) << gflops 
                  << ", \"checksum\": " << std::setprecision(4) << checksum << "}\n";
    } else {
        std::cout << "C++ Thread (" << num_threads << "T): " << elapsed_sec << "s (" << gflops << " GFLOPS)\n";
    }
    return 0;
}
