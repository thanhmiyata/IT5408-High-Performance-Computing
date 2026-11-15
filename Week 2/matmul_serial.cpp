#include "common.hpp"

int main(int argc, char* argv[]) {
    int m = 1024, n = 1024, p = 1024;
    bool json_output = false;

    if (argc >= 4) {
        m = std::atoi(argv[1]);
        n = std::atoi(argv[2]);
        p = std::atoi(argv[3]);
    }
    for (int i = 1; i < argc; ++i) {
        if (std::string(argv[i]) == "--json") json_output = true;
    }

    std::vector<double> A(static_cast<size_t>(m) * n);
    std::vector<double> B(static_cast<size_t>(n) * p);
    std::vector<double> C(static_cast<size_t>(m) * p);

    init_matrix_A(A.data(), m, n);
    init_matrix_B(B.data(), n, p);

    Timer timer;
    timer.start();
    matmul_ikj_serial(A.data(), B.data(), C.data(), m, n, p);
    double elapsed_sec = timer.stop();

    double gflops = calculate_gflops(m, n, p, elapsed_sec);
    double checksum = compute_checksum(C.data(), m, p);

    if (json_output) {
        std::cout << "{\"method\": \"serial\", \"m\": " << m << ", \"n\": " << n << ", \"p\": " << p 
                  << ", \"threads\": 1, \"time_sec\": " << std::fixed << std::setprecision(6) << elapsed_sec 
                  << ", \"gflops\": " << std::setprecision(4) << gflops 
                  << ", \"checksum\": " << std::setprecision(4) << checksum << "}\n";
    } else {
        std::cout << "Serial: " << m << "x" << n << "x" << p << " -> " << elapsed_sec << "s (" << gflops << " GFLOPS)\n";
    }
    return 0;
}
