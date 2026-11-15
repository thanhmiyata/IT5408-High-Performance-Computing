#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <climits>
#include <chrono>
#include <iomanip>
#include <algorithm>
#include <omp.h>

const int INF = 1e9;

struct TSPResult {
    int min_cost;
    std::vector<int> best_tour;
    double time_sec;
};

// Đọc dữ liệu đồ thị TSP từ file input
bool load_tsp_input(const std::string& filename, int& n, std::vector<std::vector<int>>& cost) {
    std::ifstream fin(filename);
    if (!fin.is_open()) {
        std::cerr << "Khong the mo file: " << filename << std::endl;
        return false;
    }
    fin >> n;
    cost.assign(n, std::vector<int>(n));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            fin >> cost[i][j];
        }
    }
    return true;
}

// Giải thuật Nhánh và Cận tuần tự (Branch and Bound Serial)
void tsp_branch_and_bound_serial(int n, const std::vector<std::vector<int>>& cost,
                                 int u, int count, int current_cost, int mask,
                                 std::vector<int>& path, int& best_cost, std::vector<int>& best_tour,
                                 const std::vector<int>& min_outgoing) {
    if (count == n) {
        if (cost[u][0] > 0 || u == 0) { // Quay về đỉnh xuất phát 0
            int total = current_cost + cost[u][0];
            if (total < best_cost) {
                best_cost = total;
                best_tour = path;
                best_tour.push_back(0);
            }
        }
        return;
    }

    // Tính cận dưới (Lower Bound estimation)
    int bound = current_cost;
    for (int v = 0; v < n; ++v) {
        if (!(mask & (1 << v))) {
            bound += min_outgoing[v];
        }
    }
    if (bound >= best_cost) return; // Cắt tỉa nhánh (Pruning)

    // Duyệt qua các thành phố kế tiếp
    for (int v = 0; v < n; ++v) {
        if (!(mask & (1 << v)) && cost[u][v] > 0) {
            if (current_cost + cost[u][v] < best_cost) {
                path.push_back(v);
                tsp_branch_and_bound_serial(n, cost, v, count + 1, current_cost + cost[u][v],
                                            mask | (1 << v), path, best_cost, best_tour, min_outgoing);
                path.pop_back();
            }
        }
    }
}

// Giải thuật Nhánh và Cận song song với OpenMP
TSPResult solve_tsp_parallel(int n, const std::vector<std::vector<int>>& cost, int num_threads) {
    TSPResult res;
    res.min_cost = INF;

    // Tìm cạnh nhỏ nhất đi ra từ mỗi đỉnh để tính cận dưới nhanh
    std::vector<int> min_outgoing(n, INF);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i != j && cost[i][j] < min_outgoing[i]) {
                min_outgoing[i] = cost[i][j];
            }
        }
    }

    double t_start = omp_get_wtime();

    // Nếu đồ thị rất nhỏ (n <= 4), giải tuần tự
    if (n <= 4) {
        std::vector<int> path = {0};
        tsp_branch_and_bound_serial(n, cost, 0, 1, 0, 1, path, res.min_cost, res.best_tour, min_outgoing);
        res.time_sec = omp_get_wtime() - t_start;
        return res;
    }

    // Song song hóa mức 1: Chọn đỉnh thứ hai (v1) sau đỉnh 0
    #pragma omp parallel num_threads(num_threads)
    {
        int local_best_cost = INF;
        std::vector<int> local_best_tour;

        #pragma omp for schedule(dynamic)
        for (int v1 = 1; v1 < n; ++v1) {
            // Song song hóa tiếp mức 2 nếu n >= 12
            for (int v2 = 1; v2 < n; ++v2) {
                if (v1 == v2) continue;

                std::vector<int> path = {0, v1, v2};
                int current_cost = cost[0][v1] + cost[v1][v2];
                int mask = (1 << 0) | (1 << v1) | (1 << v2);

                int branch_best = res.min_cost; // Đọc cận trên toàn cục
                tsp_branch_and_bound_serial(n, cost, v2, 3, current_cost, mask, path,
                                            branch_best, local_best_tour, min_outgoing);

                if (branch_best < res.min_cost) {
                    #pragma omp critical
                    {
                        if (branch_best < res.min_cost) {
                            res.min_cost = branch_best;
                            res.best_tour = local_best_tour;
                        }
                    }
                }
            }
        }
    }

    res.time_sec = omp_get_wtime() - t_start;
    return res;
}

int main(int argc, char* argv[]) {
    std::string filename = "tsp7.inp";
    int num_threads = 4;

    if (argc >= 2) filename = argv[1];
    if (argc >= 3) num_threads = std::atoi(argv[2]);

    int n;
    std::vector<std::vector<int>> cost;
    if (!load_tsp_input(filename, n, cost)) {
        return 1;
    }

    std::cout << "========================================================\n";
    std::cout << "  GIAI BAI TOAN NGUOI DI DU LICH (TSP - BRANCH & BOUND)\n";
    std::cout << "  File input: " << filename << " | So thanh pho (N): " << n << "\n";
    std::cout << "  So luong luong OpenMP: " << num_threads << "\n";
    std::cout << "========================================================\n";

    TSPResult result = solve_tsp_parallel(n, cost, num_threads);

    std::cout << "Chi phi nho nhat (Min Cost): " << result.min_cost << "\n";
    std::cout << "Chu trinh toi uu: ";
    for (size_t i = 0; i < result.best_tour.size(); ++i) {
        std::cout << result.best_tour[i] + 1; // 1-based indexing
        if (i + 1 < result.best_tour.size()) std::cout << " -> ";
    }
    std::cout << "\n";
    std::cout << "Thoi gian thuc thi: " << std::fixed << std::setprecision(6) 
              << result.time_sec << " giay (" << result.time_sec * 1000.0 << " ms)\n";
    std::cout << "--------------------------------------------------------\n";

    return 0;
}
