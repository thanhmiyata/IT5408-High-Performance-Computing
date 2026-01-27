# KHUNG BÁO CÁO VÀ SLIDE

## PHẦN 1: KHUNG BÁO CÁO (10 trang)

### 1. Giới thiệu (1 trang)
- **Bài toán:** Mô phỏng dòng chảy chất lưu 2D bằng Lattice Boltzmann Method.
- **Lý do chọn:**
  - Tính song song tự nhiên cao.
  - Ứng dụng thực tế: Khí động học, vi lưu, thiết kế kỹ thuật.
  - Khác biệt với FDM: Giải phương trình Boltzmann thay vì Navier-Stokes.
- **Mục tiêu:**
  - Cài đặt LBM D2Q9 tuần tự và song song MPI.
  - Đánh giá hiệu năng: MLUPS, Speedup, Efficiency.

### 2. Cơ sở lý thuyết (2 trang)
- **Phương trình Boltzmann rời rạc:**
  - Collision và Streaming.
  - Xấp xỉ BGK.
- **Mô hình D2Q9:**
  - 9 hướng vận tốc.
  - Trọng số và vector vận tốc.
- **Hàm phân bố cân bằng:**
  - Công thức Maxwell-Boltzmann rời rạc.
- **Các đại lượng vĩ mô:**
  - Mật độ ρ, vận tốc u, áp suất p.
- **Điều kiện biên:**
  - Bounce-back (tường cứng).
  - Periodic (chu kỳ).

### 3. Thuật toán số (1 trang)
- **Quy trình 2 bước:**
  1. Collision: Tính f*ᵢ từ fᵢ và fᵢᵉᑫ.
  2. Streaming: Di chuyển f*ᵢ theo hướng cᵢ.
- **Điều kiện ổn định:**
  - τ > 0.5
  - u₀ < 0.1
- **Tham số vật lý:**
  - Độ nhớt: ν = (τ - 0.5)/3
  - Số Reynolds: Re = u₀L/ν

### 4. Thiết kế song song (2 trang)
- **Chiến lược chia miền:**
  - Domain Decomposition 1D theo chiều X.
  - Mỗi tiến trình quản lý nx/P cột.
- **Ghost columns:**
  - Trao đổi 2 cột biên giữa các tiến trình lân cận.
- **Giao tiếp MPI:**
  - MPI_Irecv/Isend (non-blocking).
  - MPI_Waitall để đồng bộ.
- **Phân tích lượng giao tiếp:**
  - Dữ liệu: 2 × ny × 9 × sizeof(double) bytes mỗi bước.
  - Tỉ lệ Computation/Communication.

### 5. Cài đặt (1 trang)
- **Ngôn ngữ:** C
- **Thư viện:** MPI (OpenMPI/MPICH)
- **Cấu trúc code:**
  - `lbm_serial.c`: Bản tuần tự.
  - `lbm_mpi.c`: Bản song song MPI.
- **Các hàm chính:**
  - `KhoiTao()`: Khởi tạo fᵢ = fᵢᵉᑫ.
  - `TinhMacro()`: Tính ρ, u từ fᵢ.
  - `Collision()`: Bước va chạm.
  - `Streaming()`: Bước lan truyền.
  - `TraoDoiGhost()`: Trao đổi biên MPI.

### 6. Thực nghiệm (2 trang)
- **Môi trường:**
  - CPU: [Tên CPU, số core]
  - RAM: [Dung lượng]
  - MPI: OpenMPI [version]
- **Kịch bản test:**
  - Lưới: 256×64, 512×128, 1024×256
  - Số bước: 10000
  - Số tiến trình: 1, 2, 4, 8
- **Kết quả:**

| Lưới | P | Thời gian (s) | MLUPS | Speedup | Efficiency |
|------|---|---------------|-------|---------|------------|
| 256×64 | 1 | 1.802 | 90.94 | 1.0× | 100% |
| 256×64 | 2 | 1.365 | 120.00 | 1.32× | 66% |
| 256×64 | 4 | 0.738 | 222.13 | 2.44× | 61% |
| ... | ... | ... | ... | ... | ... |

- **Biểu đồ:**
  - Speedup vs Số tiến trình.
  - MLUPS vs Số tiến trình.
  - Efficiency vs Số tiến trình.

### 7. Phân tích và Đánh giá (1 trang)
- **Độ chính xác:**
  - So sánh checksum giữa serial và MPI.
  - Kiểm tra bảo toàn khối lượng (Σρ).
- **Hiệu năng:**
  - Speedup tốt khi lưới lớn.
  - Efficiency giảm khi tăng số tiến trình (do overhead giao tiếp).
- **Giới hạn:**
  - Khi P quá lớn so với nx, overhead áp đảo.
  - Strong scaling vs Weak scaling.

### 8. Kết luận và Hướng phát triển (0.5 trang)
- **Kết luận:**
  - Đã cài đặt thành công LBM D2Q9 tuần tự và song song.
  - Hiểu rõ cơ chế chia miền và trao đổi biên.
  - Đạt được speedup tốt trên lưới lớn.
- **Hướng phát triển:**
  - Chia miền 2D để giảm lượng giao tiếp.
  - Hybrid MPI + OpenMP.
  - Tối ưu trên GPU bằng CUDA/OpenCL.
  - Mở rộng sang 3D (D3Q19, D3Q27).

---

## PHẦN 2: KHUNG SLIDE (10-12 slide)

### Slide 1: Tiêu đề
- **Tên đề tài:** Lattice Boltzmann Method - Mô phỏng Dòng chảy 2D Song song với MPI
- **Môn học:** IT5408 - Tính toán Hiệu năng Cao
- **Sinh viên:** [Tên của bạn]
- **Giảng viên:** [Tên giảng viên]

### Slide 2: Đặt vấn đề
- **Bài toán:** Mô phỏng dòng chảy chất lưu 2D.
- **Ứng dụng:**
  - Khí động học ô tô, máy bay.
  - Thiết kế vi lưu (lab-on-a-chip).
  - Mô phỏng dòng chảy trong vật liệu xốp.
- **Thách thức:** Tính toán tốn kém → Cần song song hóa.

### Slide 3: Lattice Boltzmann Method là gì?
- **Khác biệt với FDM/FEM:**
  - FDM: Giải trực tiếp Navier-Stokes.
  - LBM: Giải phương trình Boltzmann (mô tả phân bố hạt).
- **Ưu điểm:**
  - Đơn giản: Chỉ có Collision + Streaming.
  - Song song tự nhiên: Mỗi điểm lưới độc lập.
  - Xử lý biên dễ: Bounce-back đơn giản.

### Slide 4: Mô hình D2Q9
- **Sơ đồ 9 hướng vận tốc:**
```
  6   2   5
    ↖ ↑ ↗
  3 ← 0 → 1
    ↙ ↓ ↘
  7   4   8
```
- **Trọng số:** w₀=4/9, w₁₋₄=1/9, w₅₋₈=1/36
- **Hàm phân bố:** fᵢ tại mỗi điểm, mỗi hướng.

### Slide 5: Thuật toán LBM
```
1. Collision: f*ᵢ = fᵢ - ω(fᵢ - fᵢᵉᑫ)
2. Streaming: fᵢ(x+cᵢ) = f*ᵢ(x)
3. Tính ρ, u từ fᵢ
4. Áp dụng điều kiện biên
```
- **Đơn giản:** Chỉ có phép cộng, trừ, nhân.
- **Cục bộ:** Chỉ cần dữ liệu từ điểm lân cận.

### Slide 6: Chiến lược Song song hóa
- **Chia miền 1D theo X:**
  - Mỗi tiến trình quản lý nx/P cột.
- **Ghost columns:**
  - Trao đổi 2 cột biên giữa các tiến trình.
- **Giao tiếp MPI:**
  - MPI_Irecv/Isend (non-blocking).

### Slide 7: Sơ đồ Chia miền
[Hình minh họa: Lưới 256×64 chia cho 4 tiến trình]
```
Rank 0 | Rank 1 | Rank 2 | Rank 3
64 cột | 64 cột | 64 cột | 64 cột
```
- **Trao đổi:** Rank i ↔ Rank i±1

### Slide 8: Cài đặt
- **Ngôn ngữ:** C
- **Thư viện:** MPI
- **File code:**
  - `lbm_serial.c` (tuần tự)
  - `lbm_mpi.c` (song song)
- **Tên hàm tiếng Việt không dấu:**
  - `KhoiTao`, `TinhMacro`, `Collision`, `Streaming`, `TraoDoiGhost`

### Slide 9: Kết quả Thực nghiệm
**Bảng kết quả:**
| P | Thời gian (s) | MLUPS | Speedup | Efficiency |
|---|---------------|-------|---------|------------|
| 1 | 1.802 | 90.94 | 1.0× | 100% |
| 2 | 1.365 | 120.00 | 1.32× | 66% |
| 4 | 0.738 | 222.13 | 2.44× | 61% |

**Biểu đồ Speedup:**
[Chèn biểu đồ]

### Slide 10: Phân tích Hiệu năng
- **Tại sao Speedup < P?**
  - Overhead giao tiếp MPI.
  - Overhead khởi tạo và đồng bộ.
  - Load imbalance (nếu có).
- **Khi nào hiệu quả tốt?**
  - Lưới lớn (nx >> P).
  - Tỉ lệ Computation/Communication cao.

### Slide 11: Kết luận
- **Đạt được:**
  - Cài đặt thành công LBM D2Q9 tuần tự và MPI.
  - Hiểu rõ cơ chế chia miền và trao đổi biên.
  - Speedup tốt trên lưới lớn.
- **Bài học:**
  - Tính cục bộ cao → Song song tốt.
  - Cần cân bằng giữa tính toán và giao tiếp.

### Slide 12: Hướng phát triển
- **Tối ưu hóa:**
  - Chia miền 2D.
  - Hybrid MPI + OpenMP.
  - GPU (CUDA/OpenCL).
- **Mở rộng:**
  - 3D (D3Q19, D3Q27).
  - Nhiều pha (multiphase flow).
  - Truyền nhiệt (thermal LBM).

### Slide 13: Q&A
- Cảm ơn thầy và các bạn đã lắng nghe!
- Câu hỏi thảo luận.

---

## GHI CHÚ

- **Biểu đồ cần vẽ:**
  1. Speedup vs Số tiến trình (đường lý tưởng vs thực tế).
  2. MLUPS vs Số tiến trình.
  3. Efficiency vs Số tiến trình.
  4. Sơ đồ chia miền (hình minh họa).

- **Số liệu cần đo:**
  - Chạy với nhiều kích thước lưới: 256×64, 512×128, 1024×256.
  - Chạy với nhiều số tiến trình: 1, 2, 4, 8, 16.
  - Lặp lại 3-5 lần và lấy trung bình.

- **Công cụ vẽ biểu đồ:**
  - Python (matplotlib).
  - Excel/Google Sheets.
  - Gnuplot.
