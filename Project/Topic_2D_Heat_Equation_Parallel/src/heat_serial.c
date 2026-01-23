#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

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
        exit(1);
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
    if (argc < 2) {
        printf("Usage: %s n K c\n", argv[0]);
        printf("   or: %s config.txt\n", argv[0]);
        return 1;
    }
    int n = 0;
    int K = 0;
    double c = 0.0;
    if (argc == 2) {
        read_config(argv[1], &n, &K, &c);
    } else {
        n = atoi(argv[1]);
        K = atoi(argv[2]);
        c = atof(argv[3]);
    }
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

    clock_t t0 = clock();
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
    clock_t t1 = clock();
    double elapsed = (double)(t1 - t0) / CLOCKS_PER_SEC;

    double sum = 0.0;
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            sum += at(u_old, ny, i, j);
        }
    }
    printf("Time: %.6f sec\n", elapsed);
    printf("u_mid=%.6f\n", at(u_old, ny, n / 2 + 1, n / 2 + 1));
    printf("sum=%.6f\n", sum);

    free(u_old);
    free(u_new);
    return 0;
}
