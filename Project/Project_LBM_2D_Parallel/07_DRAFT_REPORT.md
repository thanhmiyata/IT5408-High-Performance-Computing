# BÁO CÁO DỰ ÁN: MÔ PHỎNG DÒNG CHẢY 2D BẰNG LATTICE BOLTZMANN METHOD VÀ MPI

**Môn học:** IT5408 - Tính toán hiệu năng cao  
**Đề tài:** Lattice Boltzmann Method - Parallel Numerical Solution (MPI)  
**Sinh viên:** Krizpham  
**Năm:** 2026

---

## 1. GIỚI THIỆU

### 1.1. Đặt vấn đề
Mô phỏng dòng chảy chất lưu là một bài toán quan trọng trong nhiều lĩnh vực kỹ thuật như khí động học, thiết kế vi lưu, và mô phỏng vật liệu. Các phương pháp truyền thống như Finite Difference Method (FDM) và Finite Element Method (FEM) giải trực tiếp phương trình Navier-Stokes, nhưng gặp khó khăn khi xử lý biên phức tạp và song song hóa.

**Lattice Boltzmann Method (LBM)** là một phương pháp thay thế, giải phương trình Boltzmann (mô tả phân bố hạt vi mô) thay vì Navier-Stokes. LBM có những ưu điểm vượt trội:
- **Đơn giản:** Chỉ có 2 bước (Collision và Streaming).
- **Song song tự nhiên:** Mỗi điểm lưới độc lập, phù hợp với kiến trúc song song.
- **Xử lý biên dễ:** Phương pháp bounce-back đơn giản cho vật cản phức tạp.

### 1.2. Mục tiêu dự án
1. Hiểu rõ lý thuyết Lattice Boltzmann Method và mô hình D2Q9.
2. Cài đặt LBM tuần tự và song song bằng MPI.
3. Đánh giá hiệu năng: MLUPS, Speedup, Efficiency.
4. Phân tích ảnh hưởng của kích thước lưới và số tiến trình.

### 1.3. Phạm vi
- **Không gian:** 2D (hai chiều).
- **Mô hình:** D2Q9 (9 hướng vận tốc).
- **Bài toán:** Dòng chảy Poiseuille (giữa 2 tấm phẳng song song).
- **Song song hóa:** MPI với chia miền 1D theo chiều X.

---

## 2. CƠ SỞ LÝ THUYẾT

### 2.1. Phương trình Lattice Boltzmann
Phương trình LBM với xấp xỉ BGK (Bhatnagar-Gross-Krook):

```
fᵢ(x + cᵢΔt, t + Δt) = fᵢ(x, t) - ω[fᵢ(x, t) - fᵢᵉᑫ(x, t)]
```

Trong đó:
- `fᵢ(x, t)`: Hàm phân bố hạt tại vị trí x, thời điểm t, hướng i.
- `cᵢ`: Vector vận tốc rời rạc hướng i.
- `ω = 1/τ`: Tham số thư giãn.
- `fᵢᵉᑫ`: Hàm phân bố cân bằng Maxwell-Boltzmann.

### 2.2. Mô hình D2Q9
**D2Q9** có nghĩa là 2 chiều (D2) với 9 hướng vận tốc (Q9):

```
Sơ đồ hướng:
  6   2   5
    ↖ ↑ ↗
  3 ← 0 → 1
    ↙ ↓ ↘
  7   4   8
```

**Vector vận tốc:**
- c₀ = (0, 0): Đứng yên
- c₁,₂,₃,₄ = (±1, 0), (0, ±1): 4 hướng chính
- c₅,₆,₇,₈ = (±1, ±1): 4 hướng chéo

**Trọng số:**
- w₀ = 4/9, w₁₋₄ = 1/9, w₅₋₈ = 1/36

### 2.3. Hàm phân bố cân bằng
```
fᵢᵉᑫ = wᵢ × ρ × [1 + 3(cᵢ·u) + 9/2(cᵢ·u)² - 3/2(u·u)]
```

### 2.4. Các đại lượng vĩ mô
**Mật độ:**
```
ρ = Σ fᵢ
```

**Vận tốc:**
```
ρu = Σ fᵢcᵢ  →  u = (Σ fᵢcᵢ) / ρ
```

### 2.5. Điều kiện biên
- **Bounce-back (tường cứng):** fᵢ = f_opposite(i)
- **Periodic (chu kỳ):** Biên trái nối với biên phải.

---

## 3. THUẬT TOÁN VÀ CÀI ĐẶT

### 3.1. Thuật toán LBM
```
1. Khởi tạo: fᵢ = fᵢᵉᑫ với ρ=1, u=0

2. Vòng lặp thời gian (t = 0 → T):
   a. Tính ρ, u từ fᵢ
   b. Tính fᵢᵉᑫ từ ρ, u
   c. Collision: f*ᵢ = fᵢ - ω(fᵢ - fᵢᵉᑫ)
   d. Streaming: fᵢ(x+cᵢ) = f*ᵢ(x)
   e. Áp dụng điều kiện biên

3. Xuất kết quả
```

### 3.2. Thiết kế song song MPI
**Chiến lược:** Chia miền 1D theo chiều X.
- Lưới NX×NY được chia thành P khối theo chiều X.
- Mỗi tiến trình quản lý NX/P cột.

**Ghost columns:**
- Mỗi tiến trình cần 2 cột ảo (ghost_left, ghost_right) để lưu dữ liệu từ tiến trình lân cận.

**Giao tiếp MPI:**
```c
// Trao đổi biên bằng MPI_Irecv/Isend
if (rank > 0) {
  MPI_Irecv(ghost_left, ny*Q, MPI_DOUBLE, rank-1, 100, comm, &reqs[rcount++]);
  MPI_Isend(f_new + 0*ny*Q, ny*Q, MPI_DOUBLE, rank-1, 101, comm, &reqs[rcount++]);
}
if (rank < size-1) {
  MPI_Irecv(ghost_right, ny*Q, MPI_DOUBLE, rank+1, 101, comm, &reqs[rcount++]);
  MPI_Isend(f_new + (nx_local-1)*ny*Q, ny*Q, MPI_DOUBLE, rank+1, 100, comm, &reqs[rcount++]);
}
MPI_Waitall(rcount, reqs, MPI_STATUSES_IGNORE);
```

### 3.3. Cài đặt
- **Ngôn ngữ:** C
- **Thư viện:** MPI (OpenMPI)
- **File code:**
  - `lbm_serial.c`: Bản tuần tự
  - `lbm_mpi.c`: Bản song song MPI
- **Tên hàm:** Tiếng Việt không dấu (`KhoiTao`, `TinhMacro`, `Collision`, `Streaming`, `TraoDoiGhost`)

---

## 4. THỰC NGHIỆM VÀ KẾT QUẢ

### 4.1. Môi trường thực nghiệm
- **CPU:** Apple M-series / Intel Core i7 (8 cores)
- **RAM:** 16 GB
- **OS:** macOS
- **MPI:** OpenMPI 4.1.x
- **Compiler:** GCC 11.x với flag `-O2`

### 4.2. Kịch bản test
- **Lưới:** 256×64 (16,384 điểm)
- **Số bước:** 10,000
- **Tham số:** ω=1.0, u₀=0.1
- **Số tiến trình:** 1, 2, 4, 8

### 4.3. Kết quả đo lường

| Số tiến trình (P) | Thời gian (s) | MLUPS | Speedup | Efficiency |
|-------------------|---------------|-------|---------|------------|
| 1 (Serial)        | 1.802         | 90.94 | 1.00×   | 100%       |
| 2                 | 1.365         | 120.00| 1.32×   | 66%        |
| 4                 | 0.738         | 222.13| 2.44×   | 61%        |
| 8                 | 0.450         | 364.09| 4.00×   | 50%        |

**Ghi chú:**
- MLUPS = (NX × NY × NSTEPS) / (Thời gian × 10⁶)
- Speedup = T_serial / T_parallel
- Efficiency = Speedup / P × 100%

### 4.4. Biểu đồ
*(Cần vẽ biểu đồ Speedup vs P, MLUPS vs P, Efficiency vs P)*

---

## 5. PHÂN TÍCH VÀ ĐÁNH GIÁ

### 5.1. Độ chính xác
- So sánh vận tốc trung bình giữa serial và MPI: Sai số < 1e-6.
- Kiểm tra bảo toàn khối lượng: Σρ(t=0) ≈ Σρ(t=T).

### 5.2. Hiệu năng
**Ưu điểm:**
- Speedup tốt khi P nhỏ (2-4 tiến trình).
- MLUPS tăng rõ rệt khi tăng số tiến trình.

**Nhược điểm:**
- Efficiency giảm khi tăng P do overhead giao tiếp.
- Khi P=8, efficiency chỉ còn 50%.

### 5.3. Phân tích overhead
**Lượng dữ liệu trao đổi mỗi bước:**
```
Data = 2 × NY × 9 × sizeof(double) = 2 × 64 × 9 × 8 = 9,216 bytes ≈ 9 KB
```

**Tổng 10,000 bước:**
```
Total = 9 KB × 10,000 = 90 MB
```

**Tỉ lệ Computation/Communication:**
- Khi NX_local lớn → Tỉ lệ cao → Hiệu năng tốt.
- Khi NX_local nhỏ (P lớn) → Overhead áp đảo.

### 5.4. Strong Scaling
Giữ nguyên kích thước lưới, tăng số tiến trình:
- Speedup < P (không tuyến tính).
- Efficiency giảm khi tăng P.

**Nguyên nhân:**
- Overhead giao tiếp MPI.
- Overhead khởi tạo và đồng bộ.
- Định luật Amdahl (phần code tuần tự).

---

## 6. KẾT LUẬN

### 6.1. Kết quả đạt được
1. Đã cài đặt thành công LBM D2Q9 tuần tự và song song MPI.
2. Hiểu rõ cơ chế chia miền và trao đổi ghost columns.
3. Đạt được speedup 2.44× với 4 tiến trình trên lưới 256×64.
4. Đo lường được MLUPS, Speedup, Efficiency.

### 6.2. Bài học kinh nghiệm
- LBM có tính song song tự nhiên cao, phù hợp với MPI.
- Cần cân bằng giữa tính toán và giao tiếp.
- Lưới lớn cho hiệu quả song song tốt hơn.

### 6.3. Hạn chế
- Chưa tối ưu điều kiện biên vận tốc (Zou-He).
- Chưa test weak scaling.
- Chưa so sánh với chia miền 2D.

### 6.4. Hướng phát triển
1. **Tối ưu hóa:**
   - Chia miền 2D để giảm lượng giao tiếp.
   - Overlapping communication và computation.
   - Hybrid MPI + OpenMP.

2. **Mở rộng:**
   - 3D (D3Q19, D3Q27).
   - Nhiều pha (multiphase flow).
   - Truyền nhiệt (thermal LBM).
   - GPU (CUDA/OpenCL).

3. **Ứng dụng:**
   - Mô phỏng khí động học ô tô.
   - Thiết kế vi lưu (lab-on-a-chip).
   - Dòng chảy trong vật liệu xốp.

---

## 7. TÀI LIỆU THAM KHẢO

1. T. Krüger, H. Kusumaatmaja, A. Kuzmin, O. Shardt, G. Silva, and E.M. Viggen, *The Lattice Boltzmann Method: Principles and Practice*, Springer, 2017.

2. S. Chen and G.D. Doolen, "Lattice Boltzmann Method for Fluid Flows," *Annual Review of Fluid Mechanics*, vol. 30, pp. 329-364, 1998.

3. S. Succi, *The Lattice Boltzmann Equation for Fluid Dynamics and Beyond*, Oxford University Press, 2001.

4. Q. Zou and X. He, "On pressure and velocity boundary conditions for the lattice Boltzmann BGK model," *Physics of Fluids*, vol. 9, no. 6, pp. 1591-1598, 1997.

5. MPI Forum, *MPI: A Message-Passing Interface Standard*, Version 3.1, 2015. [Online]. Available: https://www.mpi-forum.org/

---

## PHỤ LỤC

### A. Code chính (lbm_serial.c)
*(Trích đoạn code quan trọng)*

```c
void Collision(double *f, double *f_new, double *rho, double *ux, double *uy, 
               int nx, int ny, double omega) {
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
```

### B. Kết quả chi tiết
*(Bảng số liệu đầy đủ từ các lần chạy thực nghiệm)*

### C. Script benchmark
*(File benchmark.sh để tự động hóa việc test)*
