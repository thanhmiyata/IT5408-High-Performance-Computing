#include <stdio.h>
#include <stdlib.h>
#include <math.h>

static double boundary_left(double y) { return 10.0; }
static double boundary_right(double y) { return 40.0; }
static double boundary_bottom(double x) { return 30.0; }
static double boundary_top(double x) { return 50.0; }
static double initial_value(double x, double y) { (void)x; (void)y; return 0.0; }

static double *alloc_grid(int nx, int ny) {
    double *a = (double *)calloc((size_t)nx * (size_t)ny, sizeof(double));
    if (!a) {
        fprintf(stderr, "Alloc failed\n");
        exit(1);
    }
    return a;
}

static inline double at(double *a, int ny, int i, int j) {
    return a[(size_t)i * (size_t)ny + (size_t)j];
}

static inline void set(double *a, int ny, int i, int j, double v) {
    a[(size_t)i * (size_t)ny + (size_t)j] = v;
}

static void init_grid(double *u, int n, double ds) {
    int nx = n + 2;
    int ny = n + 2;
    for (int i = 0; i < nx; i++) {
        double x = i * ds;
        set(u, ny, i, 0, boundary_bottom(x));
        set(u, ny, i, ny - 1, boundary_top(x));
    }
    for (int j = 0; j < ny; j++) {
        double y = j * ds;
        set(u, ny, 0, j, boundary_left(y));
        set(u, ny, nx - 1, j, boundary_right(y));
    }
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++) {
            double x = i * ds;
            double y = j * ds;
            set(u, ny, i, j, initial_value(x, y));
        }
    }
}

int main(int argc, char **argv) {
    if (argc < 4) {
        printf("Usage: %s n K c\n", argv[0]);
        return 1;
    }
    int n = atoi(argv[1]);
    int K = atoi(argv[2]);
    double c = atof(argv[3]);
    if (n <= 0 || K <= 0 || c <= 0.0) {
        fprintf(stderr, "Invalid inputs\n");
        return 1;
    }

    double ds = 1.0 / (n + 1);
    double dt = 0.25 * ds * ds / c;
    double r = c * dt / (ds * ds);

    int nx = n + 2;
    int ny = n + 2;
    double *u_old = alloc_grid(nx, ny);
    double *u_new = alloc_grid(nx, ny);

    init_grid(u_old, n, ds);
    init_grid(u_new, n, ds);

    for (int k = 0; k < K; k++) {
        for (int i = 1; i <= n; i++) {
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

    printf("u[n/2,n/2]=%.6f\n", at(u_old, ny, n / 2 + 1, n / 2 + 1));

    free(u_old);
    free(u_new);
    return 0;
}
