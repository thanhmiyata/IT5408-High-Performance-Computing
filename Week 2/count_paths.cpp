#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iomanip>
#include <omp.h>

using Matrix = std::vector<std::vector<unsigned long long>>;

// Nhân hai ma trận vuông V x V song song bằng OpenMP
Matrix matmul_parallel(const Matrix& A, const Matrix& B, int V, int num_threads) {
    Matrix C(V, std::vector<unsigned long long>(V, 0));

    // Vòng lặp IKJ tối ưu hóa Cache
    #pragma omp parallel for schedule(static) num_threads(num_threads)
    for (int i = 0; i < V; ++i) {
        for (int k = 0; k < V; ++k) {
            unsigned long long a_ik = A[i][k];
            if (a_ik == 0) continue;
            for (int j = 0; j < V; ++j) {
                C[i][j] += a_ik * B[k][j];
            }
        }
    }
    return C;
}

// Tính lũy thừa ma trận A^k bằng thuật toán lũy thừa nhị phân
Matrix matrix_power_parallel(Matrix A, int k, int V, int num_threads) {
    // Khởi tạo ma trận đơn vị I
    Matrix res(V, std::vector<unsigned long long>(V, 0));
    for (int i = 0; i < V; ++i) res[i][i] = 1;

    Matrix base = A;
    int exp = k;
    while (exp > 0) {
        if (exp % 2 == 1) {
            res = matmul_parallel(res, base, V, num_threads);
        }
        if (exp > 1) {
            base = matmul_parallel(base, base, V, num_threads);
        }
        exp /= 2;
    }
    return res;
}

// Tính tổng tất cả phần tử trong ma trận
unsigned long long sum_matrix(const Matrix& M, int V) {
    unsigned long long total = 0;
    for (int i = 0; i < V; ++i) {
        for (int j = 0; j < V; ++j) {
            total += M[i][j];
        }
    }
    return total;
}

int main(int argc, char* argv[]) {
    std::string filename = "graph01.txt";
    int k = 3;
    int num_threads = 4;
    bool json_output = false;

    if (argc >= 2) filename = argv[1];
    if (argc >= 3) k = std::atoi(argv[2]);
    if (argc >= 4) num_threads = std::atoi(argv[3]);
    for (int i = 1; i < argc; ++i) {
        if (std::string(argv[i]) == "--json") json_output = true;
    }

    std::ifstream fin(filename);
    if (!fin.is_open()) {
        std::cerr << "Khong the mo file do thi: " << filename << std::endl;
        return 1;
    }

    int V, E;
    fin >> V >> E;

    Matrix A(V, std::vector<unsigned long long>(V, 0));
    for (int i = 0; i < E; ++i) {
        int u, v;
        fin >> u >> v;
        if (u >= 1 && u <= V && v >= 1 && v <= V) {
            A[u - 1][v - 1] += 1; // 1-based to 0-based
        }
    }

    double t0 = omp_get_wtime();
    Matrix Ak = matrix_power_parallel(A, k, V, num_threads);
    unsigned long long total_paths = sum_matrix(Ak, V);
    double t1 = omp_get_wtime();
    double elapsed_sec = t1 - t0;

    if (json_output) {
        std::cout << "{\"graph\": \"" << filename << "\", \"V\": " << V << ", \"E\": " << E
                  << ", \"k\": " << k << ", \"threads\": " << num_threads 
                  << ", \"paths\": " << total_paths << ", \"time_sec\": " 
                  << std::fixed << std::setprecision(6) << elapsed_sec << "}\n";
    } else {
        std::cout << "========================================================\n";
        std::cout << "  DEM SO DUONG DI DO DAI K TRONG DO THI (OPENMP)\n";
        std::cout << "  File: " << filename << " | So dinh (V): " << V << " | So canh (E): " << E << "\n";
        std::cout << "  Do dai duong di k = " << k << " | So luong luong: " << num_threads << "\n";
        std::cout << "========================================================\n";
        std::cout << "Tong so duong di do dai " << k << ": " << total_paths << "\n";
        std::cout << "Thoi gian tinh toan: " << std::fixed << std::setprecision(6) 
                  << elapsed_sec << " giay (" << elapsed_sec * 1000.0 << " ms)\n";
        std::cout << "--------------------------------------------------------\n";
    }

    return 0;
}
