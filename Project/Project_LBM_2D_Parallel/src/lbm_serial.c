#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

// Tham số mô phỏng
#define NX  256      // Số điểm lưới theo chiều X
#define NY  64       // Số điểm lưới theo chiều Y
#define NSTEPS 30000 // Số bước thời gian (tăng để hội tụ)
#define OMEGA  1.0   // Tham số thư giãn (1/tau)
#define U0     0.1   // Vận tốc đầu vào

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
void KhoiTao(double *f, int nx, int ny) {
  int x, y, i;
  double rho = 1.0;
  double ux = 0.0;
  double uy = 0.0;
  
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
//=========================
void TinhMacro(double *f, double *rho, double *ux, double *uy, int nx, int ny) {
  int x, y, i;
  
  for (x = 0; x < nx; x++) {
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
      *(ux + x*ny + y) = vx / r;
      *(uy + x*ny + y) = vy / r;
    }
  }
}
//=========================
void Collision(double *f, double *f_new, double *rho, double *ux, double *uy, 
               int nx, int ny, double omega) {
  int x, y, i;
  
  for (x = 0; x < nx; x++) {
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
//=========================
void Streaming(double *f_new, double *f, int nx, int ny) {
  int x, y, i;
  
  for (x = 0; x < nx; x++) {
    for (y = 0; y < ny; y++) {
      for (i = 0; i < Q; i++) {
        // Tọa độ nguồn (Pull model)
        int xs = x - cx[i];
        int ys = y - cy[i];
        
        // Xử lý biên chu kỳ theo chiều X
        if (xs < 0) xs = nx - 1;
        if (xs >= nx) xs = 0;
        
        // Xử lý biên Bounce-back (vách) theo chiều Y
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
//=========================
void ApDungBienVaoCung(double *ux, int nx, int ny, double u0) {
  int y;
  // Đặt vận tốc tại biên trái
  for (y = 1; y < ny-1; y++) {
    *(ux + 0*ny + y) = u0;
    *(ux + (nx-1)*ny + y) = u0;
  }
}
//=========================
void GhiKetQua(double *rho, double *ux, double *uy, int nx, int ny, const char *filename) {
  FILE *fp = fopen(filename, "w");
  int x, y;
  
  fprintf(fp, "# x y rho ux uy\n");
  for (x = 0; x < nx; x++) {
    for (y = 0; y < ny; y++) {
      fprintf(fp, "%d %d %.6f %.6f %.6f\n", 
              x, y, 
              *(rho + x*ny + y),
              *(ux + x*ny + y),
              *(uy + x*ny + y));
    }
  }
  
  fclose(fp);
}
//=========================
int main(int argc, char **argv) {
  int nx = NX;
  int ny = NY;
  int nsteps = NSTEPS;
  double omega = OMEGA;
  double u0 = U0;
  
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
  } else if (argc > 1) {
    printf("Su dung: %s [NX] [NY] [NSTEPS] [OMEGA] [U0]\n", argv[0]);
    printf("Mac dinh: NX=%d, NY=%d, NSTEPS=%d, OMEGA=%.1f, U0=%.1f\n", NX, NY, NSTEPS, OMEGA, U0);
  }
  
  printf("=== LBM D2Q9 - Phien ban tuan tu ===\n");
  printf("Luoi: %d x %d\n", nx, ny);
  printf("So buoc: %d\n", nsteps);
  printf("Omega: %.3f\n", omega);
  printf("Van toc dau vao: %.3f\n", u0);
  
  // Cấp phát bộ nhớ
  double *f = (double *)malloc(nx * ny * Q * sizeof(double));
  double *f_new = (double *)malloc(nx * ny * Q * sizeof(double));
  double *rho = (double *)malloc(nx * ny * sizeof(double));
  double *ux = (double *)malloc(nx * ny * sizeof(double));
  double *uy = (double *)malloc(nx * ny * sizeof(double));
  
  if (!f || !f_new || !rho || !ux || !uy) {
    fprintf(stderr, "Loi cap phat bo nho!\n");
    return 1;
  }
  
  // Khởi tạo
  KhoiTao(f, nx, ny);
  
  // Bắt đầu đo thời gian
  clock_t t_start = clock();
  
  // Vòng lặp thời gian
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
  
  // Tính MLUPS
  double mlups = (double)(nx * ny * nsteps) / (elapsed * 1e6);
  
  printf("\n=== KET QUA ===\n");
  printf("Thoi gian tinh: %.3f giay\n", elapsed);
  printf("MLUPS: %.2f\n", mlups);
  
  // Tính các đại lượng vĩ mô cuối cùng
  TinhMacro(f, rho, ux, uy, nx, ny);
  
  // Tính checksum để so sánh với MPI
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
  
  // Ghi kết quả
  GhiKetQua(rho, ux, uy, nx, ny, "lbm_result.dat");
  printf("Da ghi ket qua vao: lbm_result.dat\n");
  
  // Giải phóng bộ nhớ
  free(f);
  free(f_new);
  free(rho);
  free(ux);
  free(uy);
  
  return 0;
}
