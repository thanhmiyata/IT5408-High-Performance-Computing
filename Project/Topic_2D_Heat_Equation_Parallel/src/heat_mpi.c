#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

static double boundary_left(double y) { return 10.0; }
static double boundary_right(double y) { return 40.0; }
static double boundary_bottom(double x) { return 30.0; }
static double boundary_top(double x) { return 50.0; }
static double initial_value(double x, double y) { (void)x; (void)y; return 0.0; }

static double *alloc_grid(int nx, int ny) {
    double *a = (double *)calloc((size_t)nx * (size_t)ny, sizeof(double));
    if (!a) {
        fprintf(stderr, "Alloc failed\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    return a;
}

static inline double at(double *a, int ny, int i, int j) {
    return a[(size_t)i * (size_t)ny + (size_t)j];
}

static inline void set(double *a, int ny, int i, int j, double v) {
    a[(size_t)i * (size_t)ny + (size_t)j] = v;
}

static void init_local(double *u, int n, int local_n, int start_i, double ds, int rank, int size) {
    int nx_local = local_n + 2;
    int ny = n + 2;

    for (int i = 1; i <= local_n; i++) {
        int gi = start_i + (i - 1);
        double x = gi * ds;
        set(u, ny, i, 0, boundary_bottom(x));
        set(u, ny, i, ny - 1, boundary_top(x));
    }

    if (rank == 0) {
        for (int j = 0; j < ny; j++) {
            double y = j * ds;
            set(u, ny, 0, j, boundary_left(y));
        }
    }
    if (rank == size - 1) {
        for (int j = 0; j < ny; j++) {
            double y = j * ds;
            set(u, ny, local_n + 1, j, boundary_right(y));
        }
    }

    for (int i = 1; i <= local_n; i++) {
        int gi = start_i + (i - 1);
        for (int j = 1; j <= n; j++) {
            double x = gi * ds;
            double y = j * ds;
            set(u, ny, i, j, initial_value(x, y));
        }
    }
}

static void exchange_ghosts(double *u, int local_n, int n, int rank, int size, MPI_Comm comm) {
    int ny = n + 2;
    MPI_Request reqs[4];
    int rcount = 0;

    if (rank > 0) {
        MPI_Irecv(&u[0 * ny], ny, MPI_DOUBLE, rank - 1, 100, comm, &reqs[rcount++]);
        MPI_Isend(&u[1 * ny], ny, MPI_DOUBLE, rank - 1, 101, comm, &reqs[rcount++]);
    }
    if (rank < size - 1) {
        MPI_Irecv(&u[(local_n + 1) * ny], ny, MPI_DOUBLE, rank + 1, 101, comm, &reqs[rcount++]);
        MPI_Isend(&u[local_n * ny], ny, MPI_DOUBLE, rank + 1, 100, comm, &reqs[rcount++]);
    }
    if (rcount > 0) {
        MPI_Waitall(rcount, reqs, MPI_STATUSES_IGNORE);
    }
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (argc < 4) {
        if (rank == 0) {
            printf("Usage: %s n K c\n", argv[0]);
        }
        MPI_Finalize();
        return 1;
    }

    int n = atoi(argv[1]);
    int K = atoi(argv[2]);
    double c = atof(argv[3]);
    if (n <= 0 || K <= 0 || c <= 0.0) {
        if (rank == 0) {
            fprintf(stderr, "Invalid inputs\n");
        }
        MPI_Finalize();
        return 1;
    }

    int interior_rows = n;
    int base = interior_rows / size;
    int rem = interior_rows % size;
    int local_n = base + (rank < rem ? 1 : 0);
    int start_i = 1 + rank * base + (rank < rem ? rank : rem);

    double ds = 1.0 / (n + 1);
    double dt = 0.25 * ds * ds / c;
    double r = c * dt / (ds * ds);

    int ny = n + 2;
    double *u_old = alloc_grid(local_n + 2, ny);
    double *u_new = alloc_grid(local_n + 2, ny);

    init_local(u_old, n, local_n, start_i, ds, rank, size);
    init_local(u_new, n, local_n, start_i, ds, rank, size);

    double t0 = MPI_Wtime();
    for (int k = 0; k < K; k++) {
        exchange_ghosts(u_old, local_n, n, rank, size, MPI_COMM_WORLD);
        for (int i = 1; i <= local_n; i++) {
            for (int j = 1; j <= n; j++) {
                double uij = at(u_old, ny, i, j);
                double up = at(u_old, ny, i + 1, j);
                double um = at(u_old, ny, i - 1, j);
                double vp = at(u_old, ny, i, j + 1);
                double vm = at(u_old, ny, i, j - 1);
                double next = uij + r * (up + um + vp + vm - 4.0 * uij);
                set(u_new, ny, i, j, next);
            }
        }
        double *tmp = u_old;
        u_old = u_new;
        u_new = tmp;
    }
    double t1 = MPI_Wtime();

    double local_mid = at(u_old, ny, (local_n / 2) + 1, (n / 2) + 1);
    double global_mid = 0.0;
    MPI_Reduce(&local_mid, &global_mid, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        printf("Time: %.6f sec\n", t1 - t0);
        printf("u_mid(approx)=%.6f\n", global_mid);
    }

    free(u_old);
    free(u_new);
    MPI_Finalize();
    return 0;
}
