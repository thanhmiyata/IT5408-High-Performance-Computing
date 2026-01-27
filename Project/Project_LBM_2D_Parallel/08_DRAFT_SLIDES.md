# SLIDE THUYẾT TRÌNH: LATTICE BOLTZMANN METHOD - MPI

**Thời lượng dự kiến:** 15 phút  
**Số slide:** 12-13

---

## Slide 1: TIÊU ĐỀ

**Lattice Boltzmann Method**  
**Mô phỏng Dòng chảy 2D Song song với MPI**

- **Môn học:** IT5408 - Tính toán Hiệu năng Cao
- **Sinh viên:** Krizpham
- **Giảng viên:** [Tên giảng viên]
- **Năm:** 2026

---

## Slide 2: ĐẶT VẤN ĐỀ

### Bài toán
Mô phỏng dòng chảy chất lưu 2D

### Ứng dụng thực tế
- 🚗 Khí động học ô tô, máy bay
- 🔬 Thiết kế vi lưu (lab-on-a-chip)
- 🏗️ Dòng chảy trong vật liệu xốp
- 🩺 Dòng máu trong mạch

### Thách thức
- Tính toán tốn kém → **Cần song song hóa**

---

## Slide 3: LATTICE BOLTZMANN METHOD LÀ GÌ?

### So sánh với FDM/FEM

| Phương pháp | Giải gì? | Ưu điểm | Nhược điểm |
|-------------|----------|---------|------------|
| **FDM/FEM** | Navier-Stokes trực tiếp | Chính xác cao | Phức tạp, khó song song |
| **LBM** | Phương trình Boltzmann | Đơn giản, song song tốt | Bộ nhớ lớn |

### Ưu điểm LBM
✅ **Đơn giản:** Chỉ có Collision + Streaming  
✅ **Song song tự nhiên:** Mỗi điểm lưới độc lập  
✅ **Xử lý biên dễ:** Bounce-back đơn giản  

---

## Slide 4: MÔ HÌNH D2Q9

### Sơ đồ 9 hướng vận tốc

```
  6   2   5
    ↖ ↑ ↗
  3 ← 0 → 1
    ↙ ↓ ↘
  7   4   8
```

### Thông số
- **Trọng số:** w₀ = 4/9, w₁₋₄ = 1/9, w₅₋₈ = 1/36
- **Hàm phân bố:** fᵢ tại mỗi điểm, mỗi hướng (9 giá trị)

---

## Slide 5: THUẬT TOÁN LBM

### 2 bước chính

**1. Collision (Va chạm)**
```
f*ᵢ = fᵢ - ω(fᵢ - fᵢᵉᑫ)
```
- Tính hàm phân bố sau va chạm
- ω = 1/τ: Tham số thư giãn

**2. Streaming (Lan truyền)**
```
fᵢ(x + cᵢ) = f*ᵢ(x)
```
- Di chuyển hạt theo hướng vận tốc

### Đơn giản!
- Chỉ có phép cộng, trừ, nhân
- Không có đạo hàm phức tạp

---

## Slide 6: CHIẾN LƯỢC SONG SONG HÓA

### Chia miền 1D theo X
- Lưới NX×NY chia thành P khối theo chiều X
- Mỗi tiến trình quản lý NX/P cột

### Ghost Columns
- Mỗi tiến trình cần 2 cột ảo (ghost_left, ghost_right)
- Lưu dữ liệu từ tiến trình lân cận

### Giao tiếp MPI
- `MPI_Irecv` / `MPI_Isend` (non-blocking)
- `MPI_Waitall` để đồng bộ

---

## Slide 7: SƠ ĐỒ CHIA MIỀN

### Ví dụ: Lưới 256×64 với 4 tiến trình

```
┌────────┬────────┬────────┬────────┐
│ Rank 0 │ Rank 1 │ Rank 2 │ Rank 3 │
│ 64 cột │ 64 cột │ 64 cột │ 64 cột │
└────────┴────────┴────────┴────────┘
    ↕         ↕         ↕         ↕
  Ghost    Ghost    Ghost    Ghost
```

### Trao đổi biên
- Rank i ↔ Rank i-1 (ghost_left)
- Rank i ↔ Rank i+1 (ghost_right)

---

## Slide 8: CÀI ĐẶT

### Công nghệ
- **Ngôn ngữ:** C
- **Thư viện:** MPI (OpenMPI)
- **Compiler:** GCC với flag `-O2`

### File code
- `lbm_serial.c` - Bản tuần tự
- `lbm_mpi.c` - Bản song song MPI

### Phong cách code
- Tên hàm tiếng Việt không dấu
- Truy cập mảng bằng con trỏ `*(A + i*n + j)`
- Dấu phân cách `//===` giữa các hàm

---

## Slide 9: KẾT QUẢ THỰC NGHIỆM

### Cấu hình
- **Lưới:** 256×64 (16,384 điểm)
- **Số bước:** 10,000
- **Tham số:** ω=1.0, u₀=0.1

### Kết quả

| P | Thời gian (s) | MLUPS | Speedup | Efficiency |
|---|---------------|-------|---------|------------|
| 1 | 1.802 | 90.94 | 1.0× | 100% |
| 2 | 1.365 | 120.00 | 1.32× | 66% |
| 4 | 0.738 | 222.13 | 2.44× | 61% |
| 8 | 0.450 | 364.09 | 4.00× | 50% |

**MLUPS:** Million Lattice Updates Per Second

---

## Slide 10: BIỂU ĐỒ HIỆU NĂNG

### Speedup vs Số tiến trình

```
Speedup
  8 │                    ╱ Lý tưởng
  7 │                  ╱
  6 │                ╱
  5 │              ╱
  4 │            ╱  ●  Thực tế
  3 │          ╱  ●
  2 │        ╱  ●
  1 │  ●───╱
  0 └─────────────────────
    1   2   4   6   8   P
```

*(Cần vẽ biểu đồ thực tế bằng Excel/Python)*

---

## Slide 11: PHÂN TÍCH HIỆU NĂNG

### Tại sao Speedup < P?

**Nguyên nhân:**
1. ⏱️ **Overhead giao tiếp MPI**
   - Mỗi bước: 9 KB × 2 = 18 KB
   - Tổng 10,000 bước: 180 MB

2. 🔄 **Overhead khởi tạo và đồng bộ**
   - MPI_Init, MPI_Waitall

3. 📊 **Định luật Amdahl**
   - Phần code tuần tự (đọc file, ghi kết quả)

### Khi nào hiệu quả tốt?
✅ Lưới lớn (NX >> P)  
✅ Tỉ lệ Computation/Communication cao  

---

## Slide 12: KẾT LUẬN

### Đạt được
✅ Cài đặt thành công LBM D2Q9 tuần tự và MPI  
✅ Hiểu rõ cơ chế chia miền và trao đổi ghost columns  
✅ Speedup 2.44× với 4 tiến trình  
✅ MLUPS tăng từ 90.94 → 222.13  

### Bài học
- LBM có tính song song tự nhiên cao
- Cần cân bằng giữa tính toán và giao tiếp
- Lưới lớn cho hiệu quả tốt hơn

---

## Slide 13: HƯỚNG PHÁT TRIỂN

### Tối ưu hóa
1. **Chia miền 2D** → Giảm lượng giao tiếp
2. **Hybrid MPI + OpenMP** → Tận dụng đa nhân
3. **GPU (CUDA/OpenCL)** → Tăng tốc hàng trăm lần

### Mở rộng
1. **3D** (D3Q19, D3Q27) → Mô phỏng thực tế hơn
2. **Nhiều pha** (multiphase) → Dầu-nước, khí-lỏng
3. **Truyền nhiệt** (thermal LBM) → Tản nhiệt vi mạch

### Ứng dụng
- Thiết kế khí động học ô tô
- Mô phỏng vi lưu y sinh
- Dòng chảy trong đá xốp (dầu khí)

---

## Slide 14: Q&A

### Cảm ơn thầy và các bạn đã lắng nghe!

**Câu hỏi thảo luận:**
1. Tại sao LBM phù hợp với song song hóa?
2. Ưu nhược điểm của chia miền 1D vs 2D?
3. Làm thế nào để cải thiện Efficiency?

**Liên hệ:**
- Email: [email của bạn]
- GitHub: [link repository]

---

## GHI CHÚ CHO NGƯỜI THUYẾT TRÌNH

### Thời gian phân bổ
- Slide 1-3 (Giới thiệu): 3 phút
- Slide 4-5 (Lý thuyết): 3 phút
- Slide 6-8 (Thiết kế & Cài đặt): 3 phút
- Slide 9-11 (Kết quả & Phân tích): 4 phút
- Slide 12-13 (Kết luận & Hướng phát triển): 2 phút
- Q&A: Thời gian còn lại

### Điểm nhấn
- **Slide 4:** Vẽ rõ sơ đồ D2Q9 trên bảng
- **Slide 7:** Giải thích chi tiết cách chia miền
- **Slide 9:** Nhấn mạnh MLUPS tăng rõ rệt
- **Slide 11:** Giải thích tại sao Efficiency giảm

### Demo (nếu có thời gian)
- Chạy `./lbm_serial` và `mpirun -np 4 ./lbm_mpi`
- So sánh thời gian chạy trực tiếp

### Câu hỏi dự kiến
1. **"Tại sao không dùng FDM?"**
   → LBM đơn giản hơn, song song tốt hơn, xử lý biên dễ hơn.

2. **"Chia miền 2D có tốt hơn không?"**
   → Có, giảm lượng giao tiếp (tỉ lệ diện tích/chu vi tốt hơn), nhưng code phức tạp hơn.

3. **"Tại sao Efficiency giảm khi tăng P?"**
   → Overhead giao tiếp tăng, tỉ lệ Computation/Communication giảm.

4. **"LBM có thể chạy trên GPU không?"**
   → Có, rất phù hợp với CUDA do tính song song cao và cấu trúc đều đặn.
