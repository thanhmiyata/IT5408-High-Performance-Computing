# README - Lattice Boltzmann Method (LBM) D2Q9

## 1. Mô tả
- Mô phỏng dòng chảy chất lưu 2D bằng Lattice Boltzmann Method.
- Mô hình D2Q9 (2D, 9 hướng vận tốc).
- Song song hóa bằng thư viện MPI (chia miền 1D theo chiều X).

## 2. Cấu trúc thư mục
```
Project_LBM_2D_Parallel/
├── src/
│   ├── lbm_serial.c           # Mã nguồn tuần tự
│   └── lbm_mpi.c              # Mã nguồn song song MPI
├── data/
│   ├── config_small.txt       # Cấu hình test nhỏ (tùy chọn)
│   └── config_large.txt       # Cấu hình test lớn (tùy chọn)
├── 00_OVERVIEW.md             # Tổng quan đề tài
├── 01_KIEN_THUC_CAN_HIEU.md  # Kiến thức về LBM
├── 02_MO_HINH_TOAN.md         # Mô hình toán học D2Q9
├── 03_THIET_KE_SONG_SONG_MPI.md  # Thiết kế song song
├── 04_BAO_CAO_SLIDE_OUTLINE.md   # Khung báo cáo và slide
├── 05_README.md               # File này
├── 06_TEST_DANH_GIA.md        # Hướng dẫn test và đánh giá
├── 07_DRAFT_REPORT.md         # Báo cáo chi tiết
└── 08_DRAFT_SLIDES.md         # Slide thuyết trình
```

## 3. Tham số mô phỏng

Các tham số được định nghĩa trực tiếp trong code (có thể chỉnh sửa):

```c
#define NX  256      // Số điểm lưới theo chiều X
#define NY  64       // Số điểm lưới theo chiều Y
#define NSTEPS 10000 // Số bước thời gian
#define OMEGA  1.0   // Tham số thư giãn (1/tau)
#define U0     0.1   // Vận tốc đầu vào
```

**Ý nghĩa:**
- `NX, NY`: Kích thước lưới (số điểm theo mỗi chiều).
- `NSTEPS`: Số bước thời gian mô phỏng.
- `OMEGA = 1/τ`: Tham số thư giãn (liên quan đến độ nhớt).
- `U0`: Vận tốc đặc trưng (dùng cho điều kiện biên vận tốc).

**Tính toán tham số vật lý:**
```
Độ nhớt: ν = (τ - 0.5) / 3 = (1/OMEGA - 0.5) / 3
Số Reynolds: Re = U0 × NY / ν
```

## 4. Biên dịch (Compilation)

### 4.1. Tuần tự (Serial)
```bash
gcc -O2 -o lbm_serial src/lbm_serial.c -lm
```

### 4.2. Song song (MPI)
```bash
mpicc -O2 -o lbm_mpi src/lbm_mpi.c -lm
```

**Tùy chọn tối ưu:**
- `-O2` hoặc `-O3`: Tối ưu hóa code.
- `-march=native`: Tối ưu cho CPU hiện tại.
- `-lm`: Liên kết thư viện toán học (math.h).

## 5. Chạy chương trình (Execution)

### 5.1. Tuần tự
```bash
./lbm_serial
```

**Kết quả mẫu:**
```
=== LBM D2Q9 - Phien ban tuan tu ===
Luoi: 256 x 64
So buoc: 10000
Omega: 1.000
Van toc dau vao: 0.100
Buoc 0/10000
Buoc 1000/10000
...
Buoc 9000/10000

=== KET QUA ===
Thoi gian tinh: 1.802 giay
MLUPS: 90.94
Van toc trung binh: 0.096015
Da ghi ket qua vao: lbm_result.dat
```

### 5.2. Song song (MPI)
```bash
# Chạy với 2 tiến trình
mpirun -np 2 ./lbm_mpi

# Chạy với 4 tiến trình
mpirun -np 4 ./lbm_mpi

# Chạy với 8 tiến trình
mpirun -np 8 ./lbm_mpi
```

**Kết quả mẫu:**
```
=== LBM D2Q9 - Phien ban MPI ===
So tien trinh: 4
Luoi toan cuc: 256 x 64
Luoi cuc bo (rank 0): 64 x 64
So buoc: 10000
Omega: 1.000
Buoc 0/10000
...

=== KET QUA ===
Thoi gian tinh: 0.738 giay
MLUPS: 222.13
Van toc trung binh: 0.000000
```

## 6. Chỉ số hiệu năng

### 6.1. MLUPS (Million Lattice Updates Per Second)
Số triệu điểm lưới được cập nhật mỗi giây:
```
MLUPS = (NX × NY × NSTEPS) / (Thời gian × 10^6)
```

**Ví dụ:**
- Lưới 256×64, 10000 bước, thời gian 1.802s
- MLUPS = (256 × 64 × 10000) / (1.802 × 10^6) = 90.94

### 6.2. Speedup
```
Speedup = T_serial / T_parallel
```

**Ví dụ:**
- T_serial = 1.802s
- T_parallel (4 tiến trình) = 0.738s
- Speedup = 1.802 / 0.738 = 2.44×

### 6.3. Efficiency
```
Efficiency = Speedup / P × 100%
```

**Ví dụ:**
- Speedup = 2.44×, P = 4
- Efficiency = 2.44 / 4 × 100% = 61%

## 7. File kết quả

### 7.1. lbm_result.dat (chỉ serial)
File text chứa dữ liệu đầu ra:
```
# x y rho ux uy
0 0 1.000000 0.000000 0.000000
0 1 1.000000 0.000000 0.000000
...
```

**Cột:**
- `x, y`: Tọa độ lưới
- `rho`: Mật độ
- `ux, uy`: Vận tốc theo X và Y

### 7.2. Visualize (tùy chọn)
Có thể dùng Python để vẽ trường vận tốc:
```python
import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('lbm_result.dat')
x = data[:, 0].reshape(256, 64)
y = data[:, 1].reshape(256, 64)
ux = data[:, 3].reshape(256, 64)
uy = data[:, 4].reshape(256, 64)

plt.quiver(x, y, ux, uy)
plt.xlabel('X')
plt.ylabel('Y')
plt.title('Velocity Field')
plt.show()
```

## 8. Kiểm tra kết quả

### 8.1. Bảo toàn khối lượng
Tổng mật độ phải gần như không đổi:
```
Σ ρ(t=0) ≈ Σ ρ(t=T)
```

### 8.2. So sánh Serial vs MPI
Chạy cùng tham số và so sánh:
- Vận tốc trung bình.
- Tổng mật độ.
- Giá trị tại một số điểm đại diện.

## 9. Tùy chỉnh tham số

Để thay đổi kích thước lưới hoặc số bước:
1. Mở file `src/lbm_serial.c` hoặc `src/lbm_mpi.c`.
2. Sửa các `#define` ở đầu file:
```c
#define NX  512      // Tăng lên 512
#define NY  128      // Tăng lên 128
#define NSTEPS 20000 // Tăng lên 20000 bước
```
3. Biên dịch lại.

## 10. Lưu ý

### 10.1. Điều kiện ổn định
- `OMEGA` (= 1/τ) nên trong khoảng 0.5 - 2.0.
- Nếu `OMEGA` quá lớn (τ quá nhỏ), mô phỏng sẽ không ổn định.
- Vận tốc `U0` nên < 0.1 để đảm bảo tính nén được thấp.

### 10.2. Bộ nhớ
Mỗi điểm lưới cần lưu 9 giá trị `double` (72 bytes):
- Lưới 256×64: ~1.13 MB
- Lưới 512×128: ~4.5 MB
- Lưới 1024×256: ~18 MB

### 10.3. Thời gian chạy
Phụ thuộc vào:
- Kích thước lưới (NX × NY).
- Số bước thời gian (NSTEPS).
- Hiệu năng CPU.
- Số tiến trình MPI.

**Ước lượng:**
- Lưới 256×64, 10000 bước, serial: ~2 giây.
- Lưới 512×128, 10000 bước, serial: ~8 giây.
- Lưới 1024×256, 10000 bước, serial: ~32 giây.

## 11. Troubleshooting

### 11.1. Lỗi biên dịch
```
error: mpi.h: No such file or directory
```
**Giải pháp:** Cài đặt MPI:
```bash
# macOS
brew install open-mpi

# Ubuntu/Debian
sudo apt-get install libopenmpi-dev

# CentOS/RHEL
sudo yum install openmpi-devel
```

### 11.2. Kết quả NaN hoặc Inf
**Nguyên nhân:** Mô phỏng không ổn định.
**Giải pháp:**
- Giảm `OMEGA` (tăng τ).
- Giảm `U0`.
- Giảm `NSTEPS` để kiểm tra.

### 11.3. Speedup thấp
**Nguyên nhân:**
- Lưới quá nhỏ so với số tiến trình.
- Overhead giao tiếp MPI lớn.
**Giải pháp:**
- Tăng kích thước lưới.
- Giảm số tiến trình.

## 12. Tài liệu tham khảo

1. **Sách:** T. Krüger et al., "The Lattice Boltzmann Method" (2017)
2. **Paper:** S. Chen & G.D. Doolen, "Lattice Boltzmann Method for Fluid Flows" (1998)
3. **Tutorial:** http://www.lbmethod.org/
4. **Code mẫu:** http://www.palabos.org/

## 13. Liên hệ

- **Sinh viên:** Krizpham
- **Email:** [email của bạn]
- **GitHub:** [link repository nếu có]
