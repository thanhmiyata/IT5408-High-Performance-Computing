#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

// Các giá trị mặc định cho mô phỏng
#define DEFAULT_NX  256
#define DEFAULT_NY  64
#define DEFAULT_NSTEPS 10000
#define DEFAULT_OMEGA  1.0
#define DEFAULT_U0     0.1

#define Q 9 // Mô hình D2Q9

// Trọng số cho 9 hướng vận tốc
static const double w[Q] = {
  4.0/9.0,  
  1.0/9.0,  1.0/9.0,  1.0/9.0,  1.0/9.0,  
  1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0
};

// Vector vận tốc các hướng
static const int cx[Q] = {0, 1, 0, -1,  0, 1, -1, -1,  1};
static const int cy[Q] = {0, 0, 1,  0, -1, 1,  1, -1, -1};

// Hướng đối diện (dùng cho bounce-back)
static const int opposite[Q] = {0, 3, 4, 1, 2, 7, 8, 5, 6};

/**
 * Đọc cấu hình từ file (thực hiện bởi rank 0).
 */
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

/**
 * Khởi tạo cục bộ cho từng tiến trình.
 * @param f: Mảng hàm phân phối cục bộ.
 * @param nx_local: Số cột lưới mà tiến trình này xử lý.
 * @param ny: Số dòng lưới.
 */
void KhoiTaoCucBo(double *f, int nx_local, int ny) {
  int x, y, i;
  double rho = 1.0, ux = 0.0, uy = 0.0;
  for (x = 0; x < nx_local; x++) {
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

/**
 * Tính các biến vĩ mô cục bộ trên miền của tiến trình.
 */
void TinhMacroCucBo(double *f, double *rho, double *ux, double *uy, int nx_local, int ny) {
  int x, y, i;
  for (x = 0; x < nx_local; x++) {
    for (y = 0; y < ny; y++) {
      double r = 0.0, vx = 0.0, vy = 0.0;
      for (i = 0; i < Q; i++) {
        double fi = *(f + (x*ny + y)*Q + i);
        r += fi; vx += fi * cx[i]; vy += fi * cy[i];
      }
      *(rho + x*ny + y) = r;
      if (r > 1e-10) {
        *(ux + x*ny + y) = vx / r;
        *(uy + x*ny + y) = vy / r;
      } else {
        *(ux + x*ny + y) = 0.0; *(uy + x*ny + y) = 0.0;
      }
    }
  }
}

/**
 * Bước va chạm (Collision) cục bộ.
 */
void CollisionCucBo(double *f, double *f_new, double *rho, double *ux, double *uy, int nx_local, int ny, double omega) {
  int x, y, i;
  for (x = 0; x < nx_local; x++) {
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

/**
 * Áp dụng điều kiện biên cục bộ.
 * Chỉ rank 0 (biên trái) và rank cuối (biên phải) mới áp dụng biên vận tốc.
 */
void ApDungBienCucBo(double *ux, int nx_local, int ny, double u0, int rank, int size) {
  int y;
  if (rank == 0) {
    for (y = 1; y < ny-1; y++) *(ux + 0*ny + y) = u0;
  }
  if (rank == size - 1) {
    for (y = 1; y < ny-1; y++) *(ux + (nx_local-1)*ny + y) = u0;
  }
}

/**
 * Bước lan truyền (Streaming) cục bộ.
 * Sử dụng dữ liệu từ các vùng đệm (ghost cells) nhận được từ các tiến trình lân cận.
 */
void StreamingCucBo(double *f_new, double *f, int nx_local, int ny, double *ghost_left, double *ghost_right) {
  int x, y, i;
  for (x = 0; x < nx_local; x++) {
    for (y = 0; y < ny; y++) {
      for (i = 0; i < Q; i++) {
        int xs = x - cx[i], ys = y - cy[i];
        // Biên tường (trên/dưới)
        if (ys < 0 || ys >= ny) {
          int opp = opposite[i];
          *(f + (x*ny + y)*Q + i) = *(f_new + (x*ny + y)*Q + opp);
          continue;
        }
        // Xử lý lan truyền từ tiến trình lân cận (trái/phải)
        if (xs < 0) {
          *(f + (x*ny + y)*Q + i) = *(ghost_left + y*Q + i);
        } else if (xs >= nx_local) {
          *(f + (x*ny + y)*Q + i) = *(ghost_right + y*Q + i);
        } else {
          *(f + (x*ny + y)*Q + i) = *(f_new + (xs*ny + ys)*Q + i);
        }
      }
    }
  }
}

/**
 * Trao đổi các cột biên (Ghost Columns) giữa các tiến trình lân cận.
 * Sử dụng giao tiếp không chặn (non-blocking) để tối ưu hiệu năng.
 */
void TraoDoiGhost(double *f_new, int nx_local, int ny, double *ghost_left, double *ghost_right, int rank, int size, MPI_Comm comm) {
  MPI_Request reqs[4];
  int rcount = 0;
  int left = (rank - 1 + size) % size;   // Rank bên trái (vòng lặp chu kỳ)
  int right = (rank + 1) % size;        // Rank bên phải (vòng lặp chu kỳ)

  // Nhận dữ liệu vào vùng đệm ghost
  MPI_Irecv(ghost_left, ny*Q, MPI_DOUBLE, left, 100, comm, &reqs[rcount++]);
  MPI_Irecv(ghost_right, ny*Q, MPI_DOUBLE, right, 101, comm, &reqs[rcount++]);
  
  // Gửi các cột biên thực tế cho lân cận
  MPI_Isend(f_new + 0*ny*Q, ny*Q, MPI_DOUBLE, left, 101, comm, &reqs[rcount++]);
  MPI_Isend(f_new + (nx_local-1)*ny*Q, ny*Q, MPI_DOUBLE, right, 100, comm, &reqs[rcount++]);
  
  // Đợi cho đến khi tất cả các giao tiếp hoàn tất
  MPI_Waitall(rcount, reqs, MPI_STATUSES_IGNORE);
}

int main(int argc, char **argv) {
  int rank, size;
  // Khởi tạo môi trường MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  int nx, ny, nsteps;
  double omega, u0;
  const char *config_file = "config.txt";

  // Phân tích đối số dòng lệnh
  if (argc >= 3 && strcmp(argv[1], "-c") == 0) {
    config_file = argv[2];
  } else if (argc > 1) {
    if (rank == 0) fprintf(stderr, "Loi: Cu phap khong hop le. Su dung: mpirun -np <P> %s [-c <file_config>]\n", argv[0]);
    MPI_Finalize(); return 1;
  }

  int status = 0;
  // Rank 0 đọc cấu hình và gửi cho các rank khác
  if (rank == 0) {
    printf("Dang doc cau hinh tu: %s\n", config_file);
    status = DocConfig(config_file, &nx, &ny, &nsteps, &omega, &u0);
    if (status == -1) fprintf(stderr, "Loi: Khong tim thay file config '%s'\n", config_file);
    else if (status == -2) fprintf(stderr, "Loi: File '%s' khong chua du lieu hop le\n", config_file);
  }

  // Phát tán (Broadcast) cấu hình cho mọi tiến trình
  MPI_Bcast(&status, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (status != 0) { MPI_Finalize(); return 1; }
  MPI_Bcast(&nx, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&ny, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nsteps, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&omega, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&u0, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  // Phân chia miền lưới (Domain Decomposition) theo cột
  int nx_local = nx / size;
  int remainder = nx % size;
  if (rank < remainder) nx_local++; // Xử lý nếu nx không chia hết cho size

  if (rank == 0) {
    printf("\n=== LBM D2Q9 - Phien ban MPI ===\n");
    printf("So tien trinh: %d, Luoi: %d x %d, So buoc: %d, Omega: %.3f, u0: %.3f\n", size, nx, ny, nsteps, omega, u0);
  }

  // Cấp phát bộ nhớ cục bộ
  double *f = (double *)malloc(nx_local * ny * Q * sizeof(double));
  double *f_new = (double *)malloc(nx_local * ny * Q * sizeof(double));
  double *rho = (double *)malloc(nx_local * ny * sizeof(double));
  double *ux = (double *)malloc(nx_local * ny * sizeof(double));
  double *uy = (double *)malloc(nx_local * ny * sizeof(double));
  double *ghost_left = (double *)malloc(ny * Q * sizeof(double));
  double *ghost_right = (double *)malloc(ny * Q * sizeof(double));

  if (!f || !f_new || !rho || !ux || !uy || !ghost_left || !ghost_right) {
    fprintf(stderr, "Rank %d: Loi cap phat bo nho!\n", rank);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  // Khởi tạo
  KhoiTaoCucBo(f, nx_local, ny);

  double t_start = MPI_Wtime();
  int step;

  // Vòng lặp mô phỏng song song
  for (step = 0; step < nsteps; step++) {
    TinhMacroCucBo(f, rho, ux, uy, nx_local, ny);
    ApDungBienCucBo(ux, nx_local, ny, u0, rank, size);
    CollisionCucBo(f, f_new, rho, ux, uy, nx_local, ny, omega);
    
    // Bước quan trọng: Trao đổi dữ liệu biên giữa các tiến trình lân cận
    TraoDoiGhost(f_new, nx_local, ny, ghost_left, ghost_right, rank, size, MPI_COMM_WORLD);
    
    StreamingCucBo(f_new, f, nx_local, ny, ghost_left, ghost_right);

    if (rank == 0 && step % 1000 == 0) {
      printf("Buoc %d/%d\n", step, nsteps);
    }
  }

  double t_end = MPI_Wtime();
  double elapsed = t_end - t_start;
  double mlups = (double)(nx * ny * nsteps) / (elapsed * 1e6);

  // Tính checksum để kiểm tra tính đúng đắn (so sánh với bản tuần tự)
  TinhMacroCucBo(f, rho, ux, uy, nx_local, ny);
  double local_sum_rho = 0.0, local_sum_ux = 0.0, local_sum_uy = 0.0, local_avg_ux = 0.0;
  int x, y;
  for (x = 0; x < nx_local; x++) {
    for (y = 0; y < ny; y++) {
      local_sum_rho += *(rho + x*ny + y);
      local_sum_ux += *(ux + x*ny + y);
      local_sum_uy += *(uy + x*ny + y);
      if (y > 0 && y < ny-1) {
        local_avg_ux += *(ux + x*ny + y);
      }
    }
  }

  double global_sum_rho = 0.0, global_sum_ux = 0.0, global_sum_uy = 0.0, global_avg_ux = 0.0;
  
  // Thu thập dữ liệu từ tất cả các tiến trình về rank 0
  MPI_Reduce(&local_sum_rho, &global_sum_rho, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&local_sum_ux, &global_sum_ux, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&local_sum_uy, &global_sum_uy, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&local_avg_ux, &global_avg_ux, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if (rank == 0) {
    double avg_ux = global_avg_ux / (nx * (ny-2));
    printf("\n=== KET QUA ===\n");
    printf("Thoi gian tinh: %.3f giay\n", elapsed);
    printf("MLUPS: %.2f\n", mlups);
    printf("Van toc trung binh: %.6f\n", avg_ux);
    printf("\n=== CHECKSUM (de so sanh voi Serial) ===\n");
    printf("Sum(rho): %.10f\n", global_sum_rho);
    printf("Sum(ux):  %.10f\n", global_sum_ux);
    printf("Sum(uy):  %.10f\n", global_sum_uy);
  }

  // Giải phóng bộ nhớ và kết thúc MPI
  free(f);
  free(f_new);
  free(rho);
  free(ux);
  free(uy);
  free(ghost_left);
  free(ghost_right);
  
  MPI_Finalize();
  return 0;
}
