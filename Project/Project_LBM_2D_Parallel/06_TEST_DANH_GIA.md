# HƯỚNG DẪN TEST VÀ ĐÁNH GIÁ HIỆU NĂNG

## 1. Mục tiêu thử nghiệm

### 1.1. Kiểm tra tính đúng đắn
- So sánh kết quả giữa bản tuần tự và song song.
- Kiểm tra bảo toàn khối lượng (Σρ).
- Kiểm tra tính ổn định của mô phỏng.

### 1.2. Đánh giá hiệu năng
- Đo thời gian chạy với các kích thước lưới khác nhau.
- Tính toán Speedup và Efficiency.
- Đo MLUPS (Million Lattice Updates Per Second).
- Phân tích strong scaling và weak scaling.

## 2. Kiểm tra tính đúng đắn

### 2.1. So sánh Serial vs MPI
Chạy cùng một bộ tham số và so sánh kết quả:

```bash
# Serial
./lbm_serial > serial_output.txt

# MPI với 4 tiến trình
mpirun -np 4 ./lbm_mpi > mpi_output.txt
```

**So sánh:**
- Vận tốc trung bình (avg_ux).
- Tổng mật độ (Σρ).
- Giá trị tại một số điểm đại diện.

**Kỳ vọng:**
- Sai số < 1e-6 (do làm tròn số thực).
- Nếu sai số lớn → Có lỗi trong code MPI.

### 2.2. Bảo toàn khối lượng
Tổng mật độ phải gần như không đổi:
```
Σ ρ(t=0) ≈ Σ ρ(t=T)
```

**Cách kiểm tra:**
Thêm code in ra tổng mật độ mỗi 1000 bước:
```c
if (step % 1000 == 0) {
  double total_rho = 0.0;
  for (x = 0; x < nx; x++) {
    for (y = 0; y < ny; y++) {
      total_rho += *(rho + x*ny + y);
    }
  }
  printf("Step %d: Total rho = %.6f\n", step, total_rho);
}
```

**Kỳ vọng:**
- Total rho ≈ NX × NY (vì ρ khởi tạo = 1.0).
- Biến thiên < 0.1%.

### 2.3. Kiểm tra ổn định
Nếu mô phỏng không ổn định, sẽ xuất hiện:
- NaN (Not a Number).
- Inf (Infinity).
- Giá trị tăng/giảm vô hạn.

**Nguyên nhân:**
- `OMEGA` quá lớn (τ quá nhỏ, gần 0.5).
- `U0` quá lớn (> 0.1).

**Giải pháp:**
- Giảm `OMEGA` xuống 0.8 - 1.2.
- Giảm `U0` xuống 0.05 - 0.1.

## 3. Đánh giá thời gian chạy

### 3.1. Kịch bản test
Chạy với nhiều kích thước lưới và số tiến trình:

| Lưới | Số điểm | Số bước | Mục đích |
|------|---------|---------|----------|
| 128×32 | 4,096 | 5,000 | Test nhanh |
| 256×64 | 16,384 | 10,000 | Test chuẩn |
| 512×128 | 65,536 | 10,000 | Test lớn |
| 1024×256 | 262,144 | 10,000 | Test rất lớn |

Số tiến trình: 1, 2, 4, 8, 16 (nếu có đủ core).

### 3.2. Script chạy tự động
Tạo file `benchmark.sh`:

```bash
#!/bin/bash

# Biên dịch
gcc -O2 -o lbm_serial src/lbm_serial.c -lm
mpicc -O2 -o lbm_mpi src/lbm_mpi.c -lm

# Chạy serial
echo "=== SERIAL ==="
./lbm_serial

# Chạy MPI với nhiều số tiến trình
for np in 2 4 8; do
  echo "=== MPI with $np processes ==="
  mpirun -np $np ./lbm_mpi
done
```

Chạy:
```bash
chmod +x benchmark.sh
./benchmark.sh > benchmark_results.txt
```

### 3.3. Lặp lại nhiều lần
Để kết quả chính xác, chạy mỗi cấu hình 3-5 lần và lấy trung bình:

```bash
for i in {1..5}; do
  echo "Run $i"
  mpirun -np 4 ./lbm_mpi
done
```

## 4. Tính toán chỉ số hiệu năng

### 4.1. MLUPS
```
MLUPS = (NX × NY × NSTEPS) / (Thời gian × 10^6)
```

**Ví dụ:**
- Lưới 256×64, 10000 bước
- Serial: 1.802s → MLUPS = 90.94
- MPI (4p): 0.738s → MLUPS = 222.13

### 4.2. Speedup
```
Speedup(P) = T_serial / T_parallel(P)
```

**Ví dụ:**
- T_serial = 1.802s
- T_parallel(4) = 0.738s
- Speedup(4) = 1.802 / 0.738 = 2.44×

**Speedup lý tưởng:** Speedup(P) = P
**Thực tế:** Speedup(P) < P do overhead.

### 4.3. Efficiency
```
Efficiency(P) = Speedup(P) / P × 100%
```

**Ví dụ:**
- Speedup(4) = 2.44×
- Efficiency(4) = 2.44 / 4 × 100% = 61%

**Hiệu năng tốt:** Efficiency > 70%
**Hiệu năng chấp nhận được:** Efficiency > 50%

## 5. Phân tích Strong Scaling

### 5.1. Định nghĩa
**Strong scaling:** Giữ nguyên kích thước bài toán, tăng số tiến trình.

### 5.2. Thực hiện
Chạy cùng lưới 512×128 với P = 1, 2, 4, 8:

| P | Thời gian (s) | Speedup | Efficiency |
|---|---------------|---------|------------|
| 1 | T₁ | 1.0× | 100% |
| 2 | T₂ | T₁/T₂ | (T₁/T₂)/2 × 100% |
| 4 | T₄ | T₁/T₄ | (T₁/T₄)/4 × 100% |
| 8 | T₈ | T₁/T₈ | (T₁/T₈)/8 × 100% |

### 5.3. Vẽ biểu đồ
- **Trục X:** Số tiến trình (P).
- **Trục Y:** Speedup.
- **Đường lý tưởng:** y = x (Speedup = P).
- **Đường thực tế:** Dữ liệu đo được.

## 6. Phân tích Weak Scaling

### 6.1. Định nghĩa
**Weak scaling:** Tăng kích thước bài toán tỉ lệ với số tiến trình.

### 6.2. Thực hiện
Giữ nguyên số điểm mỗi tiến trình (ví dụ: 64×64 = 4096 điểm/tiến trình):

| P | Lưới | Tổng điểm | Thời gian (s) | Efficiency |
|---|------|-----------|---------------|------------|
| 1 | 64×64 | 4,096 | T₁ | 100% |
| 2 | 128×64 | 8,192 | T₂ | T₁/T₂ × 100% |
| 4 | 256×64 | 16,384 | T₄ | T₁/T₄ × 100% |
| 8 | 512×64 | 32,768 | T₈ | T₁/T₈ × 100% |

**Kỳ vọng:** Thời gian gần như không đổi (Efficiency ≈ 100%).

## 7. Phân tích overhead

### 7.1. Thành phần thời gian
```
T_total = T_computation + T_communication + T_overhead
```

Trong đó:
- `T_computation`: Thời gian tính toán (Collision, Streaming).
- `T_communication`: Thời gian trao đổi ghost columns.
- `T_overhead`: Thời gian khởi tạo MPI, đồng bộ, etc.

### 7.2. Ước lượng lượng giao tiếp
Mỗi bước thời gian, mỗi tiến trình gửi/nhận:
```
Data = 2 × NY × 9 × sizeof(double) bytes
```

**Ví dụ:** NY=64, double=8 bytes
```
Data = 2 × 64 × 9 × 8 = 9,216 bytes ≈ 9 KB
```

Tổng 10,000 bước:
```
Total = 9 KB × 10,000 = 90 MB
```

### 7.3. Tỉ lệ Computation/Communication
```
Ratio = (NX_local × NY × 9 × Flops) / (2 × NY × 9 × Bandwidth)
```

Khi `NX_local` lớn → Ratio cao → Hiệu năng tốt.

## 8. Bảng kết quả mẫu

### 8.1. Strong Scaling (Lưới 256×64)

| P | Thời gian (s) | MLUPS | Speedup | Efficiency |
|---|---------------|-------|---------|------------|
| 1 | 1.802 | 90.94 | 1.00× | 100% |
| 2 | 1.365 | 120.00 | 1.32× | 66% |
| 4 | 0.738 | 222.13 | 2.44× | 61% |
| 8 | 0.450 | 364.09 | 4.00× | 50% |

### 8.2. Weak Scaling (4096 điểm/tiến trình)

| P | Lưới | Thời gian (s) | Efficiency |
|---|------|---------------|------------|
| 1 | 64×64 | 0.180 | 100% |
| 2 | 128×64 | 0.195 | 92% |
| 4 | 256×64 | 0.210 | 86% |
| 8 | 512×64 | 0.230 | 78% |

## 9. Kết luận và Khuyến nghị

### 9.1. Kết luận
- LBM có tính song song tự nhiên cao.
- Speedup tốt khi lưới lớn (NX >> P).
- Efficiency giảm khi tăng số tiến trình do overhead giao tiếp.

### 9.2. Khuyến nghị
- **Lưới nhỏ (< 256×64):** Dùng serial hoặc 2-4 tiến trình.
- **Lưới trung bình (512×128):** Dùng 4-8 tiến trình.
- **Lưới lớn (1024×256+):** Dùng 8-16 tiến trình.

### 9.3. Tối ưu hóa
- Sử dụng chia miền 2D để giảm lượng giao tiếp.
- Kết hợp MPI + OpenMP (hybrid).
- Tối ưu trên GPU bằng CUDA.

## 10. Checklist

- [ ] Biên dịch code với `-O2` hoặc `-O3`.
- [ ] Chạy serial và ghi lại thời gian.
- [ ] Chạy MPI với P = 2, 4, 8 và ghi lại thời gian.
- [ ] Tính Speedup và Efficiency.
- [ ] Vẽ biểu đồ Speedup vs P.
- [ ] Kiểm tra bảo toàn khối lượng.
- [ ] So sánh kết quả serial vs MPI.
- [ ] Phân tích overhead giao tiếp.
- [ ] Viết báo cáo kết quả.