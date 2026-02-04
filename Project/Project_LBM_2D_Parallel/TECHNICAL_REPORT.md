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
├── data/
│   ├── config_small.txt    # Cấu hình 128×32
│   └── config_large.txt    # Cấu hình 512×128
├── visualize_results.py    # Vẽ trường vận tốc, mật độ
├── plot_performance.py     # Vẽ biểu đồ hiệu năng
└── benchmark.sh            # Script chạy thử nghiệm
```

### 3.1. Các hàm chính trong code

| Hàm | Mô tả |
|-----|-------|
| `KhoiTao()` / `KhoiTaoCucBo()` | Khởi tạo hàm phân bố fᵢ = fᵢᵉᑫ |
| `TinhMacro()` / `TinhMacroCucBo()` | Tính ρ, uₓ, uᵧ từ fᵢ |
| `Collision()` / `CollisionCucBo()` | Bước va chạm BGK |
| `Streaming()` / `StreamingCucBo()` | Bước lan truyền (Pull model) |
| `TraoDoiGhost()` | Trao đổi ghost columns (chỉ MPI) |
| `GhiKetQua()` | Xuất file kết quả |

---

## 4. Dữ liệu Test và Tham số

### 4.1. Cấu hình thử nghiệm
| Tham số | Giá trị mặc định | Ghi chú |
|---------|------------------|---------|
| **NX × NY** | 256 × 64 | Kích thước lưới |
| **NSTEPS** | 10,000 | Số bước thời gian |
| **OMEGA (ω)** | 1.0 | Tham số thư giãn |
| **U0** | 0.1 | Vận tốc biên đầu vào |

### 4.2. Lệnh biên dịch và chạy
```bash
# Biên dịch
gcc -O2 -o lbm_serial src/lbm_serial.c -lm
mpicc -O2 -o lbm_mpi src/lbm_mpi.c -lm

# Chạy Serial
./lbm_serial 256 64 10000 1.0 0.1

# Chạy MPI với 4 tiến trình
mpirun -np 4 ./lbm_mpi 256 64 10000 1.0 0.1

# Visualize kết quả
python3 visualize_results.py
```

---

## 5. Kết quả Thực nghiệm

### 5.1. Hiệu năng tính toán (Strong Scaling)
**Cấu hình:** Lưới 256×64 = 16,384 điểm, 10,000 bước thời gian

| Số tiến trình (P) | Thời gian (s) | MLUPS | Speedup | Efficiency |
|-------------------|---------------|-------|---------|------------|
| **1 (Serial)** | 1.802 | 90.94 | 1.00× | 100% |
| **2** | 1.365 | 120.00 | 1.32× | 66% |
| **4** | 0.738 | 222.13 | 2.44× | 61% |
| **8** | 0.450 | 364.09 | 4.00× | 50% |

**Nhận xét:**
- MLUPS tăng gần **4×** từ 1 → 8 tiến trình
- Speedup tăng tuyến tính ở vùng 1-4 tiến trình
- Efficiency giảm dần do overhead giao tiếp MPI

### 5.2. Kiểm tra Tính đúng đắn

| Đại lượng | Serial | MPI | Sai lệch |
|-----------|--------|-----|----------|
| **Avg uₓ** | 0.096015 | 0.095958 | < 0.06% |

**Nhận xét:**
- Sai lệch < 0.1% → **Chấp nhận được**
- Tổng mật độ được bảo toàn tốt

### 5.3. Độ chính xác Vật lý
- **Profile vận tốc:** Dạng parabolic đúng với lý thuyết Poiseuille
- **Mật độ trung bình:** ~1.0 (bảo toàn khối lượng)
- Không xuất hiện nhiễu số hay gián đoạn tại các biên chia miền MPI

---

## 6. Kết quả Visualization

### 6.1. Trường vận tốc
- Contour plot hiển thị độ lớn vận tốc |u|
- Vector plot hiển thị hướng dòng chảy
- Vận tốc cao ở giữa kênh, thấp ở biên

### 6.2. Profile vận tốc theo Y
- So sánh LBM với lý thuyết Poiseuille
- Sai số nhỏ ở vùng biên do điều kiện bounce-back

### 6.3. Biểu đồ hiệu năng
- Speedup vs Số tiến trình
- MLUPS vs Số tiến trình  
- Efficiency vs Số tiến trình

---

## 7. Kết luận

### 7.1. Đạt được
✅ Cài đặt thành công LBM D2Q9 tuần tự và song song MPI  
✅ Hiểu rõ cơ chế chia miền và trao đổi ghost columns  
✅ Speedup **4×** với 8 tiến trình, Efficiency 50%  
✅ Kết quả vật lý khớp với lý thuyết Poiseuille  
✅ Code tự viết, không dùng thư viện LBM có sẵn

### 7.2. Hạn chế
⚠️ Efficiency giảm khi P tăng (do overhead MPI trên desktop)  
⚠️ Chưa triển khai điều kiện biên Zou-He hoàn chỉnh  
⚠️ Chưa thực hiện weak scaling

### 7.3. Hướng phát triển
- Chia miền 2D để giảm lượng giao tiếp
- Hybrid MPI + OpenMP
- Tối ưu trên GPU bằng CUDA/OpenCL
- Mở rộng sang 3D (D3Q19, D3Q27)

---

## 8. Tài liệu Tham khảo

1. T. Krüger et al., *"The Lattice Boltzmann Method: Principles and Practice"* (2017)
2. S. Chen & G.D. Doolen, *"Lattice Boltzmann Method for Fluid Flows"*, Annual Review of Fluid Mechanics (1998)
3. Q. Zou & X. He, *"On pressure and velocity boundary conditions for the lattice Boltzmann BGK model"*, Physics of Fluids (1997)
4. Tutorial: http://www.lbmethod.org/
5. Code mẫu: http://www.palabos.org/

---

*Ghi chú: Mã nguồn sử dụng mô hình Pull Streaming để đảm bảo tính nhất quán trên môi trường bộ nhớ phân tán.*
