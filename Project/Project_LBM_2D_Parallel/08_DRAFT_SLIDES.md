# GIẢI THUẬT VÀ CÀI ĐẶT CHƯƠNG TRÌNH SONG SONG
## BÀI TOÁN CÓ SỰ PHỤ THUỘC DỮ LIỆU: LATTICE BOLTZMANN METHOD (D2Q9)

**Sinh viên:** Krizpham  
**Môn học:** IT5408 - Tính toán hiệu năng cao  
**Năm:** 2026

---

## Slide 2: Mô hình toán học và Phương pháp giải

- **Phương trình Lattice Boltzmann (LBE):**
  $f_i(x + c_i \Delta t, t + \Delta t) = f_i(x, t) - \omega[f_i(x, t) - f_i^{eq}(x, t)]$
- **Công thức giải:**
  - Khởi tạo giá trị ban đầu: $f_i(x, 0) = f_i^{eq}(\rho=1, u=0)$
  - Tại bước $n+1$:
    1. **Collision:** $f^*_i(x, t) = f_i(x, t) - \omega(f_i(x, t) - f_i^{eq}(x, t))$
    2. **Streaming:** $f_i(x + c_i \Delta t, t + \Delta t) = f^*_i(x, t)$
- **Mô hình Stencil D2Q9:**
  *(Chèn ảnh lbm_d2q9_correct_stencil_diagram.png)*

---

## Slide 3: Sự phụ thuộc dữ liệu

- Trong LBM, sự phụ thuộc dữ liệu xảy ra ở bước **Streaming**:
  - Tại vị trí $x$, hàm phân bố $f_i$ ở thời điểm tiếp theo phụ thuộc vào giá trị sau va chạm $f^*_i$ từ điểm láng giềng nằm ở hướng ngược lại với $c_i$.
  - Cụ thể: $f_i(x, t+1)$ nhận giá trị từ điểm $x - c_i \Delta t$.
- **Đối với chương trình tuần tự:** Xử lý điều kiện biên (Bounce-back tại tường, Periodic tại biên X).
- **Đối với chương trình song song:** Cần trao đổi các hàm phân bố biên giữa các CPU (Ghost Columns) trước khi thực hiện Streaming.

---

## Slide 4: Các đại lượng vĩ mô và Ký hiệu

- $f_i$: Hàm phân bố hạt theo hướng $i$ ($i=0 \dots 8$).
- $c_i$: Vector vận tốc rời rạc hướng $i$.
- $w_i$: Trọng số lý thuyết ($w_0=4/9, w_{1-4}=1/9, w_{5-8}=1/36$).
- **Đại lượng vĩ mô:**
  - Mật độ: $\rho = \sum_{i=0}^{8} f_i$
  - Vận tốc: $u = \frac{1}{\rho} \sum_{i=0}^{8} f_i c_i$
- **Hàm phân bố cân bằng ($f_i^{eq}$):**
  $f_i^{eq} = w_i \rho \left(1 + 3(c_i \cdot u) + \frac{9}{2}(c_i \cdot u)^2 - \frac{3}{2}u^2\right)$

---

## Slide 5: Cài đặt thực tế: Collision & TinhMacro

```c
// Trích xuất từ lbm_mpi.c
void CollisionCucBo(double *f, double *f_new, double *rho, double *ux, double *uy, int nx_local, int ny, double omega) {
    for (int x = 0; x < nx_local; x++) {
        for (int y = 0; y < ny; y++) {
            double r = *(rho + x*ny + y);
            double usqr = ux[x*ny+y]*ux[x*ny+y] + uy[x*ny+y]*uy[x*ny+y];
            for (int i = 0; i < Q; i++) {
                double cu = cx[i]*ux[x*ny+y] + cy[i]*uy[x*ny+y];
                double feq = w[i] * r * (1.0 + 3.0*cu + 4.5*cu*cu - 1.5*usqr);
                double fi = *(f + (x*ny + y)*Q + i);
                *(f_new + (x*ny + y)*Q + i) = fi - omega * (fi - feq);
            }
        }
    }
}
```

---

## Slide 6: Cài đặt thực tế: Streaming cục bộ

```c
// Trích xuất từ lbm_mpi.c
void StreamingCucBo(double *f_new, double *f, int nx_local, int ny, double *ghost_left, double *ghost_right, int rank, int size) {
    for (int x = 0; x < nx_local; x++) {
        for (int y = 0; y < ny; y++) {
            for (int i = 0; i < Q; i++) {
                int xn = x + cx[i]; int yn = y + cy[i];
                if (yn < 0 || yn >= ny) { // Bounce-back biên Y
                    *(f + (x*ny + y)*Q + opposite[i]) += *(f_new + (x*ny + y)*Q + i);
                } else if (xn < 0) { // Lấy từ Ghost Left
                    if (rank > 0) *(f + (x*ny + yn)*Q + i) += *(ghost_left + yn*Q + i);
                    else *(f + ((nx_local-1)*ny + yn)*Q + i) += *(f_new + (x*ny + y)*Q + i); // Periodic
                } else if (xn >= nx_local) { // Lấy từ Ghost Right
                    if (rank < size-1) *(f + (x*ny + yn)*Q + i) += *(ghost_right + yn*Q + i);
                    else *(f + (0*ny + yn)*Q + i) += *(f_new + (x*ny + y)*Q + i); // Periodic
                } else {
                    *(f + (xn*ny + yn)*Q + i) += *(f_new + (x*ny + y)*Q + i);
                }
            }
        }
    }
}
```

---

## Slide 7: Giải thuật song song SPMD

- **Mô hình:** Single Program Multiple Data.
- **Chiến lược Chia miền (Domain Decomposition):**
  - Chia 1D theo cột (chiều X).
  - Mỗi tiến trình quản lý $nx\_local \approx NX/NP$ cột.
  - Vùng đệm biên: **Ghost Columns** có kích thước $NY \times Q$.
- **Giao tiếp:** Trao đổi Ghost Columns giữa các tiến trình láng giềng bằng `MPI_Isend` / `MPI_Irecv`.

---

## Slide 8: 5 Bước chính trong Giải thuật Song song

- **B1: Khởi tạo:** Cấp phát bộ nhớ ($f, f_{new}, \rho, u$, ghost buffers) và đặt trạng thái ban đầu.
- **B2: Chia miền:** Tính toán phạm vi $nx\_local$ cho mỗi Rank (xử lý phần dư $NX \% NP$).
- **B3: Phân tán:** Từng Rank tự khởi tạo vùng của mình hoặc Rank 0 phát tán dữ liệu.
- **B4: Lặp tính toán:**
  - Tính đại lượng vĩ mô $\to$ Collision $\to$ **Trao đổi biên** $\to$ Streaming.
- **B5: Tổng hợp:** Rank 0 thu thập vận tốc trung bình và ghi kết quả (Gather).

---

## Slide 9: Chi tiết Truyền thông (B4.1)

- Sử dụng cơ chế truyền thông không khóa (**Non-blocking**) để tối ưu hóa thời gian chờ.
- **Trao đổi Ghost Columns:**
  - Rank $i$ gửi cột $0$ sang Rank $i-1$, nhận vào `ghost_left`.
  - Rank $i$ gửi cột cuối sang Rank $i+1$, nhận vào `ghost_right`.
- **Đồng bộ:** Sử dụng `MPI_Waitall` trước khi thực hiện bước Streaming để đảm bảo dữ liệu vùng biên đã sẵn sàng.

---

## Slide 10: Code MPI Trao đổi biên

```c
void TraoDoiGhost(double *f_new, int nx_local, int ny, double *ghost_left, double *ghost_right, int rank, int size, MPI_Comm comm) {
    MPI_Request reqs[4]; int rcount = 0;
    if (rank > 0) { // Trao đổi với Rank trái
        MPI_Irecv(ghost_left, ny*Q, MPI_DOUBLE, rank-1, 100, comm, &reqs[rcount++]);
        MPI_Isend(f_new + 0, ny*Q, MPI_DOUBLE, rank-1, 101, comm, &reqs[rcount++]);
    }
    if (rank < size-1) { // Trao đổi với Rank phải
        MPI_Irecv(ghost_right, ny*Q, MPI_DOUBLE, rank+1, 101, comm, &reqs[rcount++]);
        MPI_Isend(f_new + (nx_local-1)*ny*Q, ny*Q, MPI_DOUBLE, rank+1, 100, comm, &reqs[rcount++]);
    }
    if (rcount > 0) MPI_Waitall(rcount, reqs, MPI_STATUSES_IGNORE);
}
```

---

## Slide 11: Kết quả thực nghiệm (Lưới 256x64)

- **Thông số:** $\omega=1.0$, $Step=10,000$, $\nu = 0.1667$.
- **Hiệu năng thực tế:**

| Số CPU (P) | Thời gian (s) | MLUPS | Speedup | Efficiency |
|---|---|---|---|---|
| 1 (Serial) | 1.978 | 82.85 | 1.00x | 100% |
| 2 | 0.993 | 165.02 | 1.99x | 99.5% |
| 4 | 0.609 | 268.88 | 3.24x | 81.0% |

- **Nhận xét:** Speedup đạt gần như lý tưởng khi P=2 và duy trì mức tốt (>80%) khi P=4.

---

## Slide 12: Đánh giá độ chính xác & Visualization

- **Sai số:** Vận tốc trung bình giữa bản MPI và Serial khớp nhau hoàn toàn (sai số máy tính).
- **Trường vận tốc:** Đạt trạng thái ổn định với profile Parabol chuẩn Poiseuille.
- *(Chèn hình ảnh velocity_field.png và velocity_profile.png)*

---

## Slide 13: Kết luận

- **Ưu điểm:** LBM D2Q9 rất phù hợp với MPI nhờ tính cục bộ của bước Collision.
- **Hạn chế:** Hiệu năng bị giới hạn bởi `Overhead` truyền thông khi kích thước lưới cục bộ quá nhỏ.
- **Hướng phát triển:** 
  - Triển khai chia miền 2D.
  - Tối ưu hóa bằng cách chồng lắp tính toán và truyền thông (overlapping).
  - Chuyển đổi sang kiến trúc GPU (CUDA).
