#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#define DEFAULT_NX  256
#define DEFAULT_NY  64
#define DEFAULT_NSTEPS 10000
#define DEFAULT_OMEGA  1.0
#define DEFAULT_U0     0.1

#define Q 9

static const double w[Q] = {
  4.0/9.0,  
  1.0/9.0,  1.0/9.0,  1.0/9.0,  1.0/9.0,  
  1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0
};

static const int cx[Q] = {0, 1, 0, -1,  0, 1, -1, -1,  1};
static const int cy[Q] = {0, 0, 1,  0, -1, 1,  1, -1, -1};

static const int opposite[Q] = {0, 3, 4, 1, 2, 7, 8, 5, 6};

int DocConfig(const char *filename, int *nx, int *ny, int *nsteps, double *omega, double *u0) {
  FILE *fp = fopen(filename, "r");
  if (!fp) return -1;
  
  char line[256];
  int found = 0;
  while (fgets(line, sizeof(line), fp)) {
    if (line[0] == '#' || line[0] == '\n' || line[0] == '\r') continue;
    char key[64], value[64];
    if (sscanf(line, "%63[^=]=%63s", key, value) == 2) {
      char *k = key;
      while (*k == ' ' || *k == '\t') k++;
      if (strcmp(k, "nx") == 0) { *nx = atoi(value); found++; }
      else if (strcmp(k, "ny") == 0) { *ny = atoi(value); found++; }
      else if (strcmp(k, "nsteps") == 0) { *nsteps = atoi(value); found++; }
      else if (strcmp(k, "omega") == 0) { *omega = atof(value); found++; }
      else if (strcmp(k, "u0") == 0) { *u0 = atof(value); found++; }
    }
  }
  fclose(fp);
  return (found > 0) ? 0 : -2;
}

void KhoiTao(double *f, int nx, int ny) {
  int x, y, i;
  double rho = 1.0, ux = 0.0, uy = 0.0;
  for (x = 0; x < nx; x++) {
    for (y = 0; y < ny; y++) {
      for (i = 0; i < Q; i++) {
        double cu = cx[i]*ux + cy[i]*uy;
        double usqr = ux*ux + uy*uy;
        double feq = w[i] * rho * (1.0 + 3.0*cu + 4.5*cu*cu - 1.5*usqr);
        *(f + (x*ny + y)*Q + i) = feq;
      }
    }
  }
}

void TinhMacro(double *f, double *rho, double *ux, double *uy, int nx, int ny) {
  int x, y, i;
  for (x = 0; x < nx; x++) {
    for (y = 0; y < ny; y++) {
      double r = 0.0, vx = 0.0, vy = 0.0;
      for (i = 0; i < Q; i++) {
        double fi = *(f + (x*ny + y)*Q + i);
        r += fi; vx += fi * cx[i]; vy += fi * cy[i];
      }
      *(rho + x*ny + y) = r;
      *(ux + x*ny + y) = vx / r;
      *(uy + x*ny + y) = vy / r;
    }
  }
}

void Collision(double *f, double *f_new, double *rho, double *ux, double *uy, int nx, int ny, double omega) {
  int x, y, i;
  for (x = 0; x < nx; x++) {
    for (y = 0; y < ny; y++) {
      double r = *(rho + x*ny + y), u_x = *(ux + x*ny + y), u_y = *(uy + x*ny + y);
      double usqr = u_x*u_x + u_y*u_y;
      for (i = 0; i < Q; i++) {
        double cu = cx[i]*u_x + cy[i]*u_y;
        double feq = w[i] * r * (1.0 + 3.0*cu + 4.5*cu*cu - 1.5*usqr);
        double fi = *(f + (x*ny + y)*Q + i);
        *(f_new + (x*ny + y)*Q + i) = fi - omega * (fi - feq);
      }
    }
  }
}

void Streaming(double *f_new, double *f, int nx, int ny) {
  int x, y, i;
  for (x = 0; x < nx; x++) {
    for (y = 0; y < ny; y++) {
      for (i = 0; i < Q; i++) {
        int xs = x - cx[i], ys = y - cy[i];
        if (xs < 0) xs = nx - 1;
        if (xs >= nx) xs = 0;
        if (ys < 0 || ys >= ny) {
          int opp = opposite[i];
          *(f + (x*ny + y)*Q + i) = *(f_new + (x*ny + y)*Q + opp);
        } else {
          *(f + (x*ny + y)*Q + i) = *(f_new + (xs*ny + ys)*Q + i);
        }
      }
    }
  }
}

void ApDungBienVaoCung(double *ux, int nx, int ny, double u0) {
  int y;
  for (y = 1; y < ny-1; y++) {
    *(ux + 0*ny + y) = u0;
    *(ux + (nx-1)*ny + y) = u0;
  }
}

void GhiKetQua(double *rho, double *ux, double *uy, int nx, int ny, const char *filename) {
  FILE *fp = fopen(filename, "w");
  int x, y;
  fprintf(fp, "# x y rho ux uy\n");
  for (x = 0; x < nx; x++) {
    for (y = 0; y < ny; y++) {
      fprintf(fp, "%d %d %.6f %.6f %.6f\n", x, y, *(rho + x*ny + y), *(ux + x*ny + y), *(uy + x*ny + y));
    }
  }
  fclose(fp);
}

int main(int argc, char **argv) {
  int nx, ny, nsteps;
  double omega, u0;
  const char *config_file = "config.txt";
  if (argc >= 3 && strcmp(argv[1], "-c") == 0) {
    config_file = argv[2];
  } else if (argc > 1) {
    fprintf(stderr, "Loi: Cu phap khong hop le. Su dung: %s [-c <file_cau_hinh>]\n", argv[0]);
    return 1;
  }
  int status = DocConfig(config_file, &nx, &ny, &nsteps, &omega, &u0);
  if (status == -1) {
    fprintf(stderr, "Loi: Khong tim thay file config '%s'\n", config_file);
    return 1;
  } else if (status == -2) {
    fprintf(stderr, "Loi: File '%s' khong chua du lieu cau hinh hop le (nx, ny, nsteps, omega, u0)\n", config_file);
    return 1;
  }
  printf("\n=== LBM D2Q9 - Code tuan tu ===\n");
  printf("Luoi: %d x %d\n", nx, ny);
  printf("So buoc: %d\n", nsteps);
  printf("Omega: %.3f\n", omega);
  printf("Van toc dau vao: %.3f\n", u0);
  double *f = (double *)malloc(nx * ny * Q * sizeof(double));
  double *f_new = (double *)malloc(nx * ny * Q * sizeof(double));
  double *rho = (double *)malloc(nx * ny * sizeof(double));
  double *ux = (double *)malloc(nx * ny * sizeof(double));
  double *uy = (double *)malloc(nx * ny * sizeof(double));
  if (!f || !f_new || !rho || !ux || !uy) {
    fprintf(stderr, "Loi cap phat bo nho!\n");
    return 1;
  }
  KhoiTao(f, nx, ny);
  clock_t t_start = clock();
  int step;
  for (step = 0; step < nsteps; step++) {
    TinhMacro(f, rho, ux, uy, nx, ny);
    ApDungBienVaoCung(ux, nx, ny, u0);
    Collision(f, f_new, rho, ux, uy, nx, ny, omega);
    Streaming(f_new, f, nx, ny);
    if (step % 1000 == 0) {
      printf("Buoc %d/%d\n", step, nsteps);
    }
  }
  clock_t t_end = clock();
  double elapsed = (double)(t_end - t_start) / CLOCKS_PER_SEC;
  double mlups = (double)(nx * ny * nsteps) / (elapsed * 1e6);
  printf("\n=== KET QUA ===\n");
  printf("Thoi gian tinh: %.3f giay\n", elapsed);
  printf("MLUPS: %.2f\n", mlups);
  TinhMacro(f, rho, ux, uy, nx, ny);
  double sum_rho = 0.0;
  double sum_ux = 0.0;
  double sum_uy = 0.0;
  double avg_ux = 0.0;
  int x, y;
  for (x = 0; x < nx; x++) {
    for (y = 0; y < ny; y++) {
      sum_rho += *(rho + x*ny + y);
      sum_ux += *(ux + x*ny + y);
      sum_uy += *(uy + x*ny + y);
      if (y > 0 && y < ny-1) {
        avg_ux += *(ux + x*ny + y);
      }
    }
  }
  avg_ux /= (nx * (ny-2));
  printf("Van toc trung binh: %.6f\n", avg_ux);
  printf("\n=== CHECKSUM (de so sanh voi MPI) ===\n");
  printf("Sum(rho): %.10f\n", sum_rho);
  printf("Sum(ux):  %.10f\n", sum_ux);
  printf("Sum(uy):  %.10f\n", sum_uy);
  GhiKetQua(rho, ux, uy, nx, ny, "lbm_result.dat");
  printf("Da ghi ket qua vao: lbm_result.dat\n");
  free(f);
  free(f_new);
  free(rho);
  free(ux);
  free(uy);
  return 0;
}
