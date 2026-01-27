#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

// Tham số mặc định
static double g_left = 10.0;
static double g_right = 40.0;
static double g_bottom = 30.0;
static double g_top = 50.0;
static double g_init = 0.0;
//=========================
static double DieuKienBienTrai(double y) { 
  (void)y; 
  return g_left; 
}
//=========================
static double DieuKienBienPhai(double y) { 
  (void)y; 
  return g_right; 
}
//=========================
static double DieuKienBienDuoi(double x) { 
  (void)x; 
  return g_bottom; 
}
//=========================
static double DieuKienBienTren(double x) { 
  (void)x; 
  return g_top; 
}
//=========================
static double GiaTriKhoiTao(double x, double y) { 
  (void)x; 
  (void)y; 
  return g_init; 
}
//=========================
static void XoaKyTuXuongDong(char *s) {
  size_t len = strlen(s);
  if (len > 0 && (s[len - 1] == '\n' || s[len - 1] == '\r')) {
    s[len - 1] = '\0';
  }
}
//=========================
static void DocFileCauHinh(const char *path, int *n, int *K, double *c) {
  FILE *f = fopen(path, "r");
  if (!f) {
    fprintf(stderr, "Khong the mo file cau hinh: %s\n", path);
    exit(1);
  }
  char line[256];
  while (fgets(line, sizeof(line), f)) {
    XoaKyTuXuongDong(line);
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
//=========================
static double *CapPhatLuoi(int nx, int ny) {
  double *a = (double *)calloc((size_t)nx * (size_t)ny, sizeof(double));
  if (!a) {
    fprintf(stderr, "Cap phat bo nho that bai\n");
    exit(1);
  }
  return a;
}
//=========================
static void KhoiTaoLuoi(double *u, int n, double ds) {
  int i, j;
  int nx = n + 2;
  int ny = n + 2;
  
  // Điều kiện biên dưới và trên
  for (i = 0; i < nx; i++) {
    double x = i * ds;
    *(u + i*ny + 0) = DieuKienBienDuoi(x);
    *(u + i*ny + (ny-1)) = DieuKienBienTren(x);
  }
  
  // Điều kiện biên trái và phải
  for (j = 0; j < ny; j++) {
    double y = j * ds;
    *(u + 0*ny + j) = DieuKienBienTrai(y);
    *(u + (nx-1)*ny + j) = DieuKienBienPhai(y);
  }
  
  // Giá trị khởi tạo bên trong
  for (i = 1; i <= n; i++) {
    for (j = 1; j <= n; j++) {
      double x = i * ds;
      double y = j * ds;
      *(u + i*ny + j) = GiaTriKhoiTao(x, y);
    }
  }
}
//=========================
static void TinhSaiPhan(double *u_old, double *u_new, int n, double r) {
  int i, j;
  int ny = n + 2;
  
  for (i = 1; i <= n; i++) {
    for (j = 1; j <= n; j++) {
      double uij = *(u_old + i*ny + j);
      double up  = *(u_old + (i+1)*ny + j);
      double um  = *(u_old + (i-1)*ny + j);
      double vp  = *(u_old + i*ny + (j+1));
      double vm  = *(u_old + i*ny + (j-1));
      *(u_new + i*ny + j) = uij + r * (up + um + vp + vm - 4.0 * uij);
    }
  }
}
//=========================
int main(int argc, char **argv) {
  int i, j, k;
  
  if (argc < 2) {
    printf("Cach dung: %s <file_cau_hinh.txt>\n", argv[0]);
    printf("       or: %s <n> <K> <c>\n", argv[0]);
    return 1;
  }
  
  int n = 0;
  int K = 0;
  double c = 0.0;
  
  if (argc == 2) {
    DocFileCauHinh(argv[1], &n, &K, &c);
  } else {
    n = atoi(argv[1]);
    K = atoi(argv[2]);
    c = atof(argv[3]);
  }
  
  if (n <= 0 || K <= 0 || c <= 0.0) {
    fprintf(stderr, "Tham so dau vao khong hop le\n");
    return 1;
  }
  
  double ds = 1.0 / (n + 1);
  double dt = 0.25 * ds * ds / c;
  double r = c * dt / (ds * ds);
  
  int nx = n + 2;
  int ny = n + 2;
  double *u_old = CapPhatLuoi(nx, ny);
  double *u_new = CapPhatLuoi(nx, ny);
  
  KhoiTaoLuoi(u_old, n, ds);
  KhoiTaoLuoi(u_new, n, ds);
  
  // Bắt đầu đo thời gian
  clock_t t0 = clock();
  
  // Vòng lặp thời gian
  for (k = 0; k < K; k++) {
    TinhSaiPhan(u_old, u_new, n, r);
    // Hoán đổi con trỏ
    double *tmp = u_old;
    u_old = u_new;
    u_new = tmp;
  }
  
  clock_t t1 = clock();
  double elapsed = (double)(t1 - t0) / CLOCKS_PER_SEC;
  
  // Tính tổng và giá trị giữa
  double sum = 0.0;
  for (i = 0; i < nx; i++) {
    for (j = 0; j < ny; j++) {
      sum += *(u_old + i*ny + j);
    }
  }
  
  printf("Thoi gian tinh: %.6f giay\n", elapsed);
  printf("u_mid = %.6f\n", *(u_old + (n/2+1)*ny + (n/2+1)));
  printf("sum = %.6f\n", sum);
  
  free(u_old);
  free(u_new);
  return 0;
}
