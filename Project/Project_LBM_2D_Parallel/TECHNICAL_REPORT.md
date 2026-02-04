# BÁO CÁO KỸ THUẬT: MÔ PHỎNG DÒNG CHẢY 2D BẰNG LBM VÀ MPI

**Môn học:** Tính toán hiệu năng cao (IT5408)  
**Sinh viên:** Krizpham  
**Ngày cập nhật:** 04/02/2026

---

## 1. Tổng quan Đề tài

### 1.1. Mục tiêu
Mô phỏng dòng chảy chất lưu 2D (dòng Poiseuille) bằng phương pháp **Lattice Boltzmann Method (LBM)** và song song hóa với **MPI** để đánh giá hiệu năng tính toán.

### 1.2. Tại sao chọn LBM?
| Tiêu chí | Lý do |
|----------|-------|
| **Tính song song tự nhiên** | Mỗi điểm lưới chỉ tương tác với điểm lân cận |
| **Xử lý biên dễ dàng** | Phương pháp bounce-back đơn giản |
| **Ứng dụng thực tế** | Khí động học, vi lưu, y sinh |
| **Khác biệt với FDM** | Giải phương trình Boltzmann thay vì Navier-Stokes |

---

## 2. Mô tả Kỹ thuật và Cài đặt

### 2.1. Kỹ thuật mô phỏng
- **Phương pháp:** Lattice Boltzmann Method (LBM) với mô hình **D2Q9**
- **Bài toán:** Mô phỏng dòng chảy Poiseuille trong kênh 2 chiều
- **Điều kiện biên:**
  - Trên/Dưới: Bounce-back (tường cứng)
  - Trái/Phải: Periodic (chu kỳ) + Vận tốc cưỡng ép

### 2.2. Chiến lược song song (MPI)
- **Mô hình:** SPMD (Single Program Multiple Data)
- **Chia miền:** Domain Decomposition 1D theo chiều X (chia theo cột)
- **Ghost Columns:** Trao đổi 2 cột biên giữa các tiến trình lân cận
- **Giao tiếp:** Non-blocking (`MPI_Isend`, `MPI_Irecv`, `MPI_Waitall`)

### 2.3. Hai bước tính toán chính
1. **Collision (Va chạm):** `f*ᵢ = fᵢ - ω(fᵢ - fᵢᵉᑫ)`
2. **Streaming (Lan truyền):** `fᵢ(x + cᵢ, t+1) = f*ᵢ(x, t)`

### 2.4. Cài đặt hệ thống
| Thành phần | Chi tiết |
|------------|----------|
| **Ngôn ngữ** | C |
| **Thư viện MPI** | OpenMPI |
| **Môi trường** | macOS / Linux |
| **Trình biên dịch** | `gcc` (serial), `mpicc` (MPI) |
| **Tối ưu hóa** | `-O2` |

---

## 3. Cấu trúc Mã nguồn

```
Project_LBM_2D_Parallel/
├── src/
│   ├── lbm_serial.c        # Phiên bản tuần tự
│   └── lbm_mpi.c           # Phiên bản MPI
├── config.txt              # File config mặc định
├── visualize_results.py    # Vẽ trường vận tốc, mật độ
├── plot_performance.py     # Vẽ biểu đồ hiệu năng
└── benchmark.sh            # Script chạy thử nghiệm
```

### 3.1. Các hàm chính trong code

| Hàm | Mô tả |
|-----|-------|
| `DocConfig()` | Đọc tham số từ file config |
| `KhoiTao()` / `KhoiTaoCucBo()` | Khởi tạo hàm phân bố fᵢ = fᵢᵉᑫ |
| `TinhMacro()` / `TinhMacroCucBo()` | Tính ρ, uₓ, uᵧ từ fᵢ |
| `Collision()` / `CollisionCucBo()` | Bước va chạm BGK |
| `Streaming()` / `StreamingCucBo()` | Bước lan truyền (Pull model) |
| `TraoDoiGhost()` | Trao đổi ghost columns (chỉ MPI) |
| `GhiKetQua()` | Xuất file kết quả |

---

## 4. Dữ liệu Test và Tham số

### 4.1. File Config (`config.txt`)
```
# Cấu hình LBM mặc định
nx=256
ny=64
nsteps=10000
omega=1.0
u0=0.1
```

### 4.2. Cấu hình thử nghiệm
| Tham số | Giá trị mặc định | Ghi chú |
|---------|------------------|---------|
| **NX × NY** | 256 × 64 | Kích thước lưới (16,384 điểm) |
| **NSTEPS** | 10,000 | Số bước thời gian |
| **OMEGA (ω)** | 1.0 | Tham số thư giãn |
| **U0** | 0.1 | Vận tốc biên đầu vào |

### 4.3. Lệnh biên dịch và chạy
```bash
# Biên dịch
gcc -O2 -o lbm_serial src/lbm_serial.c -lm
mpicc -O2 -o lbm_mpi src/lbm_mpi.c -lm

# Chạy Serial (tự động đọc config.txt)
./lbm_serial

# Chạy MPI với 4 tiến trình (tự động đọc config.txt)
mpirun -np 4 ./lbm_mpi

# Chạy với tham số dòng lệnh tùy chỉnh
mpirun -np 4 ./lbm_mpi 256 64 10000 1.0 0.1

# Visualize kết quả
python3 visualize_results.py
python3 plot_performance.py
```

---

## 5. Kết quả Thực nghiệm

### 5.1. Hiệu năng tính toán (Strong Scaling)
**Cấu hình:** Lưới 256×64 = 16,384 điểm, 10,000 bước thời gian

| Số tiến trình (P) | Thời gian (s) | MLUPS | Speedup | Efficiency |
|-------------------|---------------|-------|---------|------------|
| **1 (Serial)** | 1.591 | 102.95 | 1.00× | 100% |
| **2** | 0.990 | 165.45 | 1.61× | 80.5% |
| **4** | 0.621 | 263.98 | 2.56× | 64.1% |

**Công thức:**
- **MLUPS** = (NX × NY × NSTEPS) / (Time × 10⁶)
- **Speedup** = T_serial / T_parallel
- **Efficiency** = Speedup / P × 100%

**Nhận xét:**
- MLUPS tăng **2.5×** từ Serial → 4 tiến trình
- Speedup đạt **2.56×** với 4 tiến trình
- Efficiency **64%** (chấp nhận được cho desktop)

### 5.2. Kiểm tra Tính đúng đắn (Checksum)

| Đại lượng | Serial | MPI (4P) | Sai lệch |
|-----------|--------|----------|----------|
| **Sum(ρ)** | 16383.9999 | 16381.4419 | 0.016% ✅ |
| **Sum(uₓ)** | 1527.0476 | 1526.0273 | 0.067% ✅ |
| **Sum(uᵧ)** | 0.0000 | -0.0000 | 0% ✅ |
| **Avg uₓ** | 0.096015 | 0.095958 | 0.06% ✅ |

**Nhận xét:**
- Sai lệch < 0.1% → **Kết quả đáng tin cậy**
- Tổng mật độ được bảo toàn tốt
- MPI implementation chính xác

### 5.3. Độ chính xác Vật lý
- **Profile vận tốc:** Dạng parabolic đúng với lý thuyết Poiseuille
  - u_max (LBM) = 0.1403
  - u_max (Fit Parabolic) = 0.1392
  - RMSE = 0.000648
- **Mật độ trung bình:** ~1.0 (bảo toàn khối lượng)
- Không xuất hiện nhiễu số hay gián đoạn tại các biên chia miền MPI

---

## 6. Kết quả Visualization

### 6.1. Trường vận tốc (`velocity_field.png`)
- Contour plot hiển thị độ lớn vận tốc |u|
- Vector plot hiển thị hướng dòng chảy
- Vận tốc cao ở giữa kênh (màu đỏ), thấp ở biên (màu xanh)

### 6.2. Profile vận tốc (`velocity_profile.png`)
- So sánh LBM với Fit Parabolic (Poiseuille)
- Hiển thị sai số giữa LBM và đường cong lý thuyết
- Max Error = 0.001319

### 6.3. Trường mật độ (`density_field.png`)
- Mật độ dao động nhẹ quanh 1.0 (range: 0.95 - 1.07)
- Bảo toàn khối lượng tốt

### 6.4. Biểu đồ hiệu năng (`performance_charts.png`)
- Speedup vs Số tiến trình
- MLUPS vs Số tiến trình  
- Efficiency vs Số tiến trình

---

## 7. Kết luận

### 7.1. Đạt được
✅ Cài đặt thành công LBM D2Q9 tuần tự và song song MPI  
✅ Hiểu rõ cơ chế chia miền và trao đổi ghost columns  
✅ Speedup **2.56×** với 4 tiến trình, Efficiency 64%  
✅ Kết quả vật lý khớp với lý thuyết Poiseuille  
✅ Code tự viết, không dùng thư viện LBM có sẵn  
✅ Đọc tham số từ file config linh hoạt

### 7.2. Hạn chế
⚠️ Efficiency giảm khi P tăng (do overhead MPI trên desktop)  
⚠️ Chưa triển khai điều kiện biên Zou-He hoàn chỉnh  
⚠️ Chưa thực hiện weak scaling  
⚠️ Chia miền 1D chưa tối ưu cho mọi trường hợp

### 7.3. Hướng phát triển
- Chia miền **2D** để giảm lượng giao tiếp
- **Hybrid MPI + OpenMP** cho node SMP
- Tối ưu trên **GPU** bằng CUDA/OpenCL
- Mở rộng sang **3D** (D3Q19, D3Q27)

---

## 8. Tài liệu Tham khảo

1. T. Krüger et al., *"The Lattice Boltzmann Method: Principles and Practice"* (2017)
2. S. Chen & G.D. Doolen, *"Lattice Boltzmann Method for Fluid Flows"*, Annual Review of Fluid Mechanics (1998)
3. Q. Zou & X. He, *"On pressure and velocity boundary conditions for the lattice Boltzmann BGK model"*, Physics of Fluids (1997)
4. OpenMPI Documentation: https://www.open-mpi.org/doc/
5. Tutorial LBM: http://www.lbmethod.org/

---

*Ghi chú: Mã nguồn sử dụng mô hình Pull Streaming để đảm bảo tính nhất quán trên môi trường bộ nhớ phân tán.*
