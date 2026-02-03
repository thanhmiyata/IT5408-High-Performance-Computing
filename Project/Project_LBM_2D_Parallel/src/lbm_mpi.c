#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Hằng số D2Q9
#define Q 9

// Trọng số cho 9 hướng
static const double w[Q] = {
  4.0/9.0,  // 0: center
  1.0/9.0,  1.0/9.0,  1.0/9.0,  1.0/9.0,  // 1-4: cardinal
  1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0  // 5-8: diagonal
};

// Vector vận tốc rời rạc (cx, cy)
static const int cx[Q] = {0, 1, 0, -1,  0, 1, -1, -1,  1};
static const int cy[Q] = {0, 0, 1,  0, -1, 1,  1, -1, -1};

// Hướng đối diện (cho bounce-back)
static const int opposite[Q] = {0, 3, 4, 1, 2, 7, 8, 5, 6};
//=========================
void KhoiTaoCucBo(double *f, int nx_local, int ny) {
  int x, y, i;
  double rho = 1.0;
  double ux = 0.0;
  double uy = 0.0;
  
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
//=========================
void TinhMacroCucBo(double *f, double *rho, double *ux, double *uy, 
                     int nx_local, int ny) {
  int x, y, i;
  
  for (x = 0; x < nx_local; x++) {
    for (y = 0; y < ny; y++) {
      double r = 0.0;
      double vx = 0.0;
      double vy = 0.0;
      
      for (i = 0; i < Q; i++) {
        double fi = *(f + (x*ny + y)*Q + i);
        r += fi;
        vx += fi * cx[i];
        vy += fi * cy[i];
      }
      
      *(rho + x*ny + y) = r;
      if (r > 1e-10) {
        *(ux + x*ny + y) = vx / r;
        *(uy + x*ny + y) = vy / r;
      } else {
        *(ux + x*ny + y) = 0.0;
        *(uy + x*ny + y) = 0.0;
      }
    }
  }
}
//=========================
void CollisionCucBo(double *f, double *f_new, double *rho, double *ux, double *uy,
                     int nx_local, int ny, double omega) {
  int x, y, i;
  
  for (x = 0; x < nx_local; x++) {
    for (y = 0; y < ny; y++) {
      double r = *(rho + x*ny + y);
      double u_x = *(ux + x*ny + y);
      double u_y = *(uy + x*ny + y);
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
//=========================
void ApDungBienCucBo(double *ux, int nx_local, int ny, double u0, int rank, int size) {
  int y;
  // Nếu là rank đầu (chứa biên trái)
  if (rank == 0) {
    for (y = 1; y < ny-1; y++) {
      *(ux + 0*ny + y) = u0;
    }
  }
  // Nếu là rank cuối (chứa biên phải)
  if (rank == size - 1) {
    for (y = 1; y < ny-1; y++) {
      *(ux + (nx_local-1)*ny + y) = u0;
    }
  }
}
//=========================
void StreamingCucBo(double *f_new, double *f, int nx_local, int ny,
                     double *ghost_left, double *ghost_right) {
  int x, y, i;
  
  for (x = 0; x < nx_local; x++) {
    for (y = 0; y < ny; y++) {
      for (i = 0; i < Q; i++) {
        // source coordinate
        int xs = x - cx[i];
        int ys = y - cy[i];
        
        // Xử lý biên Y (bounce-back)
        if (ys < 0 || ys >= ny) {
          int opp = opposite[i];
          *(f + (x*ny + y)*Q + i) = *(f_new + (x*ny + y)*Q + opp);
          continue;
        }
        
        // Xử lý biên X (Pull from ghost cells if needed)
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
void TraoDoiGhost(double *f_new, int nx_local, int ny, 
                   double *ghost_left, double *ghost_right,
                   int rank, int size, MPI_Comm comm) {
  MPI_Request reqs[4];
  int rcount = 0;
  
  int left = (rank - 1 + size) % size;
  int right = (rank + 1) % size;

  // Nhận từ láng giềng
  MPI_Irecv(ghost_left, ny*Q, MPI_DOUBLE, left, 100, comm, &reqs[rcount++]);
  MPI_Irecv(ghost_right, ny*Q, MPI_DOUBLE, right, 101, comm, &reqs[rcount++]);
  
  // Gửi cho láng giềng
  MPI_Isend(f_new + 0*ny*Q, ny*Q, MPI_DOUBLE, left, 101, comm, &reqs[rcount++]);
  MPI_Isend(f_new + (nx_local-1)*ny*Q, ny*Q, MPI_DOUBLE, right, 100, comm, &reqs[rcount++]);
  
  MPI_Waitall(rcount, reqs, MPI_STATUSES_IGNORE);
}
//=========================
int main(int argc, char **argv) {
  int rank, size;
  
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  
  // Tham số
  int nx = 256;
  int ny = 64;
  int nsteps = 30000;  // Tăng để hội tụ tốt hơn
  double omega = 1.0;
  double u0 = 0.1;
  
  if (argc >= 6) {
    nx = atoi(argv[1]);
    ny = atoi(argv[2]);
    nsteps = atoi(argv[3]);
    omega = atof(argv[4]);
    u0 = atof(argv[5]);
  } else if (argc >= 4) {
    nx = atoi(argv[1]);
    ny = atoi(argv[2]);
    nsteps = atoi(argv[3]);
  }
  
  // Chia miền theo X
  int nx_local = nx / size;
  int remainder = nx % size;
  if (rank < remainder) nx_local++;
  
  if (rank == 0) {
    printf("=== LBM D2Q9 - Phien ban MPI ===\n");
    printf("So tien trinh: %d\n", size);
    printf("Luoi toan cuc: %d x %d\n", nx, ny);
    printf("Luoi cuc bo (rank 0): %d x %d\n", nx_local, ny);
    printf("So buoc: %d\n", nsteps);
    printf("Omega: %.3f\n", omega);
  }
  
  // Cấp phát bộ nhớ
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
  
  // Bắt đầu đo thời gian
  double t_start = MPI_Wtime();
  
  // Vòng lặp thời gian
  int step;
  for (step = 0; step < nsteps; step++) {
    TinhMacroCucBo(f, rho, ux, uy, nx_local, ny);
    ApDungBienCucBo(ux, nx_local, ny, u0, rank, size);
    CollisionCucBo(f, f_new, rho, ux, uy, nx_local, ny, omega);
    TraoDoiGhost(f_new, nx_local, ny, ghost_left, ghost_right, rank, size, MPI_COMM_WORLD);
    StreamingCucBo(f_new, f, nx_local, ny, ghost_left, ghost_right);
    
    if (rank == 0 && step % 1000 == 0) {
      printf("Buoc %d/%d\n", step, nsteps);
    }
  }
  
  double t_end = MPI_Wtime();
  double elapsed = t_end - t_start;
  
  // Tính MLUPS
  double mlups = (double)(nx * ny * nsteps) / (elapsed * 1e6);
  
  // Tính các đại lượng vĩ mô cuối cùng
  TinhMacroCucBo(f, rho, ux, uy, nx_local, ny);
  
  // Tính checksum cục bộ
  double local_sum_rho = 0.0;
  double local_sum_ux = 0.0;
  double local_sum_uy = 0.0;
  double local_avg_ux = 0.0;
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
  
  // Tổng hợp checksum từ tất cả các tiến trình
  double global_sum_rho = 0.0;
  double global_sum_ux = 0.0;
  double global_sum_uy = 0.0;
  double global_avg_ux = 0.0;
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
  
  // Giải phóng bộ nhớ
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
