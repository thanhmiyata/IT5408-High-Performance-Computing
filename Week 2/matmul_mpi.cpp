#include "common.hpp"
#include <mpi.h>

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);
    int rank, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

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

    std::vector<int> sendcounts_A(nprocs, 0), displs_A(nprocs, 0);
    std::vector<int> recvcounts_C(nprocs, 0), displs_C(nprocs, 0);
    std::vector<int> rows_per_proc(nprocs, 0);

    int base_rows = m / nprocs;
    int rem_rows = m % nprocs;
    int offset_row = 0;

    for (int r = 0; r < nprocs; ++r) {
        int r_count = base_rows + (r < rem_rows ? 1 : 0);
        rows_per_proc[r] = r_count;
        sendcounts_A[r] = r_count * n;
        displs_A[r] = offset_row * n;
        recvcounts_C[r] = r_count * p;
        displs_C[r] = offset_row * p;
        offset_row += r_count;
    }

    int local_m = rows_per_proc[rank];
    std::vector<double> A, C;
    if (rank == 0) {
        A.resize(static_cast<size_t>(m) * n);
        C.resize(static_cast<size_t>(m) * p);
        init_matrix_A(A.data(), m, n);
    }

    std::vector<double> B(static_cast<size_t>(n) * p);
    if (rank == 0) init_matrix_B(B.data(), n, p);

    std::vector<double> local_A(static_cast<size_t>(local_m) * n);
    std::vector<double> local_C(static_cast<size_t>(local_m) * p, 0.0);

    MPI_Barrier(MPI_COMM_WORLD);
    double t_start = MPI_Wtime();

    MPI_Bcast(B.data(), n * p, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatterv(A.data(), sendcounts_A.data(), displs_A.data(), MPI_DOUBLE,
                 local_A.data(), local_m * n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    for (int i = 0; i < local_m; ++i) {
        for (int k = 0; k < n; ++k) {
            double a_ik = local_A[i * n + k];
            for (int j = 0; j < p; ++j) {
                local_C[i * p + j] += a_ik * B[k * p + j];
            }
        }
    }

    MPI_Gatherv(local_C.data(), local_m * p, MPI_DOUBLE,
                C.data(), recvcounts_C.data(), displs_C.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);
    double elapsed_sec = MPI_Wtime() - t_start;
    double max_time = 0.0;
    MPI_Reduce(&elapsed_sec, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        double gflops = calculate_gflops(m, n, p, max_time);
        double checksum = compute_checksum(C.data(), m, p);
        if (json_output) {
            std::cout << "{\"method\": \"mpi\", \"m\": " << m << ", \"n\": " << n << ", \"p\": " << p 
                      << ", \"threads\": " << nprocs << ", \"time_sec\": " << std::fixed << std::setprecision(6) << max_time 
                      << ", \"gflops\": " << std::setprecision(4) << gflops 
                      << ", \"checksum\": " << std::setprecision(4) << checksum << "}\n";
        } else {
            std::cout << "MPI (" << nprocs << "P): " << max_time << "s (" << gflops << " GFLOPS)\n";
        }
    }

    MPI_Finalize();
    return 0;
}
