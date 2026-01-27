#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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
    MPI_Abort(MPI_COMM_WORLD, 1);
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
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  return a;
}
//=========================
static void KhoiTaoLuoiCucBo(double *u, int n, int local_n, int start_i, 
                              double ds, int rank, int size) {
  int i, j;
  int nx_local = local_n + 2;
  int ny = n + 2;
  
  // Điều kiện biên dưới và trên cho các hàng cục bộ
  for (i = 1; i <= local_n; i++) {
    int gi = start_i + (i - 1);
    double x = gi * ds;
    *(u + i*ny + 0) = DieuKienBienDuoi(x);
    *(u + i*ny + (ny-1)) = DieuKienBienTren(x);
  }
  
  // Điều kiện biên trái (chỉ rank 0)
  if (rank == 0) {
    for (j = 0; j < ny; j++) {
      double y = j * ds;
      *(u + 0*ny + j) = DieuKienBienTrai(y);
    }
  }
  
  // Điều kiện biên phải (chỉ rank cuối)
  if (rank == size - 1) {
    for (j = 0; j < ny; j++) {
      double y = j * ds;
      *(u + (local_n+1)*ny + j) = DieuKienBienPhai(y);
    }
  }
  
  // Giá trị khởi tạo bên trong
  for (i = 1; i <= local_n; i++) {
    int gi = start_i + (i - 1);
    for (j = 1; j <= n; j++) {
      double x = gi * ds;
      double y = j * ds;
      *(u + i*ny + j) = GiaTriKhoiTao(x, y);
    }
  }
}
//=========================
static void TraoDoiBien(double *u, int local_n, int n, int rank, int size) {
  int ny = n + 2;
  MPI_Request reqs[4];
  int rcount = 0;
  
  // Trao đổi với rank trên (rank - 1)
  if (rank > 0) {
    MPI_Irecv(u + 0*ny, ny, MPI_DOUBLE, rank - 1, 100, 
              MPI_COMM_WORLD, &reqs[rcount++]);
    MPI_Isend(u + 1*ny, ny, MPI_DOUBLE, rank - 1, 101, 
              MPI_COMM_WORLD, &reqs[rcount++]);
  }
  
  // Trao đổi với rank dưới (rank + 1)
  if (rank < size - 1) {
    MPI_Irecv(u + (local_n+1)*ny, ny, MPI_DOUBLE, rank + 1, 101, 
              MPI_COMM_WORLD, &reqs[rcount++]);
    MPI_Isend(u + local_n*ny, ny, MPI_DOUBLE, rank + 1, 100, 
              MPI_COMM_WORLD, &reqs[rcount++]);
  }
  
  if (rcount > 0) {
    MPI_Waitall(rcount, reqs, MPI_STATUSES_IGNORE);
  }
}
//=========================
static void TinhSaiPhan(double *u_old, double *u_new, int local_n, int n, double r) {
  int i, j;
  int ny = n + 2;
  
  for (i = 1; i <= local_n; i++) {
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
  int rank, size;
  
  // Khởi tạo MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  
  if (argc < 2) {
    if (rank == 0) {
      printf("Cach dung: %s <file_cau_hinh.txt>\n", argv[0]);
      printf("       or: %s <n> <K> <c>\n", argv[0]);
    }
    MPI_Finalize();
    return 1;
  }
  
  int n = 0;
  int K = 0;
  double c = 0.0;
  
  // Rank 0 đọc file cấu hình
  if (rank == 0) {
    if (argc == 2) {
      DocFileCauHinh(argv[1], &n, &K, &c);
    } else {
      n = atoi(argv[1]);
      K = atoi(argv[2]);
      c = atof(argv[3]);
    }
  }
  
  // Broadcast tham số cho tất cả các rank
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
      fprintf(stderr, "Tham so dau vao khong hop le\n");
    }
    MPI_Finalize();
    return 1;
  }
  
  // Chia miền theo hàng
  int interior_rows = n;
  int base = interior_rows / size;
  int rem = interior_rows % size;
  int local_n = base + (rank < rem ? 1 : 0);
  int start_i = 1 + rank * base + (rank < rem ? rank : rem);
  
  double ds = 1.0 / (n + 1);
  double dt = 0.25 * ds * ds / c;
  double r = c * dt / (ds * ds);
  
  int ny = n + 2;
  double *u_old = CapPhatLuoi(local_n + 2, ny);
  double *u_new = CapPhatLuoi(local_n + 2, ny);
  
  KhoiTaoLuoiCucBo(u_old, n, local_n, start_i, ds, rank, size);
  KhoiTaoLuoiCucBo(u_new, n, local_n, start_i, ds, rank, size);
  
  // Bắt đầu đo thời gian
  double t0 = MPI_Wtime();
  
  // Vòng lặp thời gian
  for (k = 0; k < K; k++) {
    TraoDoiBien(u_old, local_n, n, rank, size);
    TinhSaiPhan(u_old, u_new, local_n, n, r);
    // Hoán đổi con trỏ
    double *tmp = u_old;
    u_old = u_new;
    u_new = tmp;
  }
  
  double t1 = MPI_Wtime();
  
  // Tính tổng cục bộ
  double local_sum = 0.0;
  for (i = 0; i < local_n + 2; i++) {
    for (j = 0; j < ny; j++) {
      local_sum += *(u_old + i*ny + j);
    }
  }
  
  // Reduce tổng về rank 0
  double global_sum = 0.0;
  MPI_Reduce(&local_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  
  // Tính giá trị giữa
  double mid_val = 0.0;
  int mid_i = (n / 2) + 1;
  if (mid_i >= start_i && mid_i < start_i + local_n) {
    int li = (mid_i - start_i) + 1;
    mid_val = *(u_old + li*ny + (n/2+1));
  }
  
  double global_mid = 0.0;
  MPI_Reduce(&mid_val, &global_mid, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  
  // Rank 0 in kết quả
  if (rank == 0) {
    printf("Thoi gian tinh: %.6f giay\n", t1 - t0);
    printf("u_mid = %.6f\n", global_mid);
    printf("sum = %.6f\n", global_sum);
  }
  
  free(u_old);
  free(u_new);
  MPI_Finalize();
  return 0;
}
