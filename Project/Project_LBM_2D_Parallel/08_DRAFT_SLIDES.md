# GIẢI THUẬT VÀ CÀI ĐẶT CHƯƠNG TRÌNH SONG SONG
## BÀI TOÁN CÓ SỰ PHỤ THUỘC DỮ LIỆU: LATTICE BOLTZMANN METHOD (D2Q9)

**Sinh viên:** Krizpham  
**Môn học:** IT5408 - Tính toán hiệu năng cao  
**Năm:** 2026

---

## Slide 1: Mô hình Toán học và Phương pháp Giải

### Phương trình Lattice Boltzmann (LBE)
$$f_i(x + c_i \Delta t, t + \Delta t) = f_i(x, t) - \omega[f_i(x, t) - f_i^{eq}(x, t)]$$

### Mô hình D2Q9 - 9 hướng vận tốc
```
    6   2   5
      \ | /
    3 - 0 - 1
      / | \
    7   4   8
```

### Hai bước tính toán chính
| Bước | Công thức | Mô tả |
|------|-----------|-------|
| **Collision** | $f^*_i = f_i - \omega(f_i - f_i^{eq})$ | Va chạm (cục bộ) |
| **Streaming** | $f_i(x + c_i, t+1) = f^*_i(x, t)$ | Lan truyền (láng giềng) |

*(Chèn ảnh: lbm_d2q9_stencil.png)*

---

## Slide 2: Sự Phụ Thuộc Dữ Liệu

### Vấn đề
Trong bước **Streaming**, $f_i(x, t+1)$ phụ thuộc vào $f^*_i(x - c_i, t)$ từ **điểm láng giềng**.

### Giải pháp
| Mô hình | Xử lý biên |
|---------|------------|
| **Tuần tự** | Bounce-back (tường Y), Periodic (biên X) |
| **Song song** | Trao đổi **Ghost Columns** giữa các tiến trình |

### Minh họa Ghost Columns
```
  Rank 0          Rank 1          Rank 2
┌─────────┐     ┌─────────┐     ┌─────────┐
│ ▓▓ ... ◀│────▶│▓▓ ... ▓▓│◀───▶│... ▓▓  │
│ Ghost   │     │ Ghost   │     │ Ghost  │
└─────────┘     └─────────┘     └─────────┘
     ↑               ↑               ↑
  nx_local        nx_local        nx_local
```

---

## Slide 3: Đại Lượng Vĩ Mô và Ký Hiệu

### Ký hiệu chính
| Ký hiệu | Mô tả | Giá trị |
|---------|-------|---------|
| $f_i$ | Hàm phân bố hạt hướng $i$ | $i = 0...8$ |
| $c_i$ | Vector vận tốc rời rạc | $(c_x, c_y)$ |
| $w_i$ | Trọng số | $w_0=4/9$, $w_{1-4}=1/9$, $w_{5-8}=1/36$ |
| $\omega$ | Tham số thư giãn | $\omega = 1.0$ |

### Đại lượng vĩ mô (tính từ $f_i$)
$$\rho = \sum_{i=0}^{8} f_i \qquad \mathbf{u} = \frac{1}{\rho} \sum_{i=0}^{8} f_i \mathbf{c}_i$$

### Hàm phân bố cân bằng
$$f_i^{eq} = w_i \rho \left(1 + 3(\mathbf{c}_i \cdot \mathbf{u}) + \frac{9}{2}(\mathbf{c}_i \cdot \mathbf{u})^2 - \frac{3}{2}|\mathbf{u}|^2\right)$$

---

## Slide 4: Cài Đặt Thực Tế

### Code Collision (BGK)
```c
void CollisionCucBo(double *f, double *f_new, double *rho, 
                    double *ux, double *uy, int nx_local, 
                    int ny, double omega) {
    for (int x = 0; x < nx_local; x++) {
        for (int y = 0; y < ny; y++) {
            double r = rho[x*ny + y];
            double usqr = ux[...]*ux[...] + uy[...]*uy[...];
            for (int i = 0; i < Q; i++) {
                double cu = cx[i]*ux[...] + cy[i]*uy[...];
                double feq = w[i] * r * (1.0 + 3.0*cu 
                           + 4.5*cu*cu - 1.5*usqr);
                f_new[...] = f[...] - omega*(f[...] - feq);
            }
        }
    }
}
```
**Đặc điểm:** Hoàn toàn cục bộ - không cần giao tiếp MPI!

---

## Slide 5: Giải Thuật Song Song SPMD

### Chiến lược: Domain Decomposition 1D
```
┌───────────────────────────────────────────────┐
│  Rank 0  │  Rank 1  │  Rank 2  │  Rank 3     │
│  64 cột  │  64 cột  │  64 cột  │  64 cột     │
└───────────────────────────────────────────────┘
           ↕         ↕         ↕
      Ghost Exchange via MPI_Isend/Irecv
```

### 5 Bước chính
| Bước | Mô tả | MPI Operations |
|------|-------|----------------|
| **B1** | Khởi tạo MPI | `MPI_Init`, `MPI_Comm_rank` |
| **B2** | Cấp phát bộ nhớ cục bộ | f, rho, u, ghost_left/right |
| **B3** | Khởi tạo f = f_eq | ρ=1, u=0 |
| **B4** | **Vòng lặp:** Collision → Ghost → Streaming | `MPI_Isend/Irecv` |
| **B5** | Tổng hợp kết quả | `MPI_Reduce` |

---

## Slide 6: Chi Tiết Truyền Thông

### Code Trao đổi Ghost Columns (Non-blocking)
```c
void TraoDoiGhost(...) {
    MPI_Request reqs[4];
    int left = (rank - 1 + size) % size;
    int right = (rank + 1) % size;
    
    // Non-blocking receive
    MPI_Irecv(ghost_left, ny*Q, MPI_DOUBLE, left, 100, ...);
    MPI_Irecv(ghost_right, ny*Q, MPI_DOUBLE, right, 101, ...);
    
    // Non-blocking send
    MPI_Isend(f_new + 0, ny*Q, MPI_DOUBLE, left, 101, ...);
    MPI_Isend(f_new + (nx_local-1)*ny*Q, ny*Q, MPI_DOUBLE, right, 100, ...);
    
    MPI_Waitall(4, reqs, MPI_STATUSES_IGNORE);
}
```

### Lý do dùng Non-blocking
- Giảm thời gian chờ (overlap)
- Tránh deadlock

---

## Slide 7: Kết Quả Thực Nghiệm

### Cấu hình: Lưới 256×64, 10,000 bước, ω=1.0

| P | Thời gian (s) | MLUPS | Speedup | Efficiency |
|---|---------------|-------|---------|------------|
| **1** | 1.594 | 102.79 | 1.00× | 100% |
| **2** | 0.950 | 172.55 | **1.68×** | **84.0%** |
| **4** | 0.691 | 236.96 | **2.31×** | **57.7%** |
| **8** | 1.111 | 147.48 | **1.44×** | **17.9%** |

*(Chèn ảnh: performance_charts.png)*

### Kiểm tra đúng đắn (Checksum)
| Đại lượng | Serial | MPI (8P) | Sai lệch |
|-----------|--------|----------|----------|
| Sum(ρ) | 16383.99 | 16380.95 | 0.018% ✅ |
| Avg(uₓ) | 0.096015 | 0.096032 | 0.017% ✅ |

*(Chèn ảnh: velocity_field.png, velocity_profile.png)*

---

## Slide 8: Kết Luận

### ✅ Đã đạt được
- Cài đặt LBM D2Q9 tuần tự và song song MPI
- Speedup **2.56×** với 4 tiến trình
- Kết quả khớp lý thuyết Poiseuille (sai số < 0.1%)
- Đọc config từ file linh hoạt

### ⚠️ Hạn chế
- Efficiency giảm khi P > 4 (overhead MPI)
- Chia miền 1D chưa tối ưu cho mọi trường hợp

### 🚀 Hướng phát triển
- Chia miền **2D** để giảm lượng giao tiếp
- **Hybrid MPI + OpenMP**
- Chuyển sang **GPU (CUDA/OpenCL)**

---

## Tài liệu Tham khảo

1. T. Krüger et al., *"The Lattice Boltzmann Method"* (2017)
2. S. Chen & G.D. Doolen, *"LBM for Fluid Flows"*, Ann. Rev. Fluid Mech. (1998)
3. OpenMPI Documentation: https://www.open-mpi.org/

---

## Q&A - Cảm ơn đã lắng nghe!

**Source code:** `Project_LBM_2D_Parallel/`

| File | Mô tả |
|------|-------|
| `src/lbm_serial.c` | Phiên bản tuần tự |
| `src/lbm_mpi.c` | Phiên bản MPI |
| `config.txt` | File cấu hình |
