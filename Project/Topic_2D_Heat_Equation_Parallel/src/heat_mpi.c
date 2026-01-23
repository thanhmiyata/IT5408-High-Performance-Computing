#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static double g_left = 10.0;
static double g_right = 40.0;
static double g_bottom = 30.0;
static double g_top = 50.0;
static double g_init = 0.0;

static double boundary_left(double y) { (void)y; return g_left; }
static double boundary_right(double y) { (void)y; return g_right; }
static double boundary_bottom(double x) { (void)x; return g_bottom; }
static double boundary_top(double x) { (void)x; return g_top; }
static double initial_value(double x, double y) { (void)x; (void)y; return g_init; }

static void trim_newline(char *s) {
    size_t len = strlen(s);
    if (len > 0 && (s[len - 1] == '\n' || s[len - 1] == '\r')) {
        s[len - 1] = '\0';
    }
}

static void read_config(const char *path, int *n, int *K, double *c) {
    FILE *f = fopen(path, "r");
    if (!f) {
        fprintf(stderr, "Cannot open config: %s\n", path);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    char line[256];
    while (fgets(line, sizeof(line), f)) {
        trim_newline(line);
        if (line[0] == '#' || line[0] == '\0') {
            continue;
        }
        char key[64];
        char val[64];
        if (sscanf(line, "%63[^=]=%63s", key, val) == 2) {
            if (strcmp(key, "n") == 0) *n = atoi(val);
            else if (strcmp(key, "K") == 0) *K = atoi(val);
            else if (strcmp(key, "c") == 0) *c = atof(val);
            else if (strcmp(key, "left") == 0) g_left = atof(val);
            else if (strcmp(key, "right") == 0) g_right = atof(val);
            else if (strcmp(key, "bottom") == 0) g_bottom = atof(val);
            else if (strcmp(key, "top") == 0) g_top = atof(val);
            else if (strcmp(key, "init") == 0) g_init = atof(val);
        }
    }
    fclose(f);
}

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

    if (argc < 2) {
        if (rank == 0) {
            printf("Usage: %s n K c\n", argv[0]);
            printf("   or: %s config.txt\n", argv[0]);
        }
        MPI_Finalize();
        return 1;
    }

    int n = 0;
    int K = 0;
    double c = 0.0;
    if (rank == 0) {
        if (argc == 2) {
            read_config(argv[1], &n, &K, &c);
        } else {
            n = atoi(argv[1]);
            K = atoi(argv[2]);
            c = atof(argv[3]);
        }
    }
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&K, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&c, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&g_left, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&g_right, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&g_bottom, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&g_top, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&g_init, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
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

    double local_sum = 0.0;
    for (int i = 0; i < local_n + 2; i++) {
        for (int j = 0; j < ny; j++) {
            local_sum += at(u_old, ny, i, j);
        }
    }
    double global_sum = 0.0;
    MPI_Reduce(&local_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    double mid_val = 0.0;
    int mid_i = (n / 2) + 1;
    if (mid_i >= start_i && mid_i < start_i + local_n) {
        int li = (mid_i - start_i) + 1;
        mid_val = at(u_old, ny, li, (n / 2) + 1);
    }
    double global_mid = 0.0;
    MPI_Reduce(&mid_val, &global_mid, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        printf("Time: %.6f sec\n", t1 - t0);
        printf("u_mid=%.6f\n", global_mid);
        printf("sum=%.6f\n", global_sum);
    }

    free(u_old);
    free(u_new);
    MPI_Finalize();
    return 0;
}
