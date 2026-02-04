# BÁO CÁO KỸ THUẬT: MÔ PHỎNG DÒNG CHẢY 2D BẰNG LBM VÀ MPI

**Môn học:** Tính toán hiệu năng cao (IT5408)  
**Sinh viên:** Krizpham  
**Ngày cập nhật:** 04/02/2026

## 1. Mô tả Kỹ thuật và Cài đặt

### Kỹ thuật mô phỏng
- **Phương pháp:** Lattice Boltzmann Method (LBM) với mô hình **D2Q9**.
- **Bài toán:** Mô phỏng dòng chảy Poiseuille trong kênh 2 chiều (ổn định).
- **Cơ chế song song:** Chia miền một chiều (1D Domain Decomposition) theo cột (trục X). Sử dụng các cột **Ghost Columns** để trao đổi dữ liệu biên giữa các tiến trình lân cận qua giao tiếp không khóa (**Non-blocking communication**: `MPI_Isend`, `MPI_Irecv`).
- **Nghiệm lý thuyết:** Profile vận tốc có dạng Parabol.

### Cài đặt hệ thống
- **Ngôn ngữ:** C (chuẩn C99).
- **Thư viện:** OpenMPI.
- **Biên dịch:** Sử dụng `gcc` cho phiên bản tuần tự và `mpicc` cho phiên bản song song với cờ tối ưu hóa `-O2`.

## 2. Source Code

### Cấu trúc thư mục
```text
Project_LBM_2D_Parallel/
├── src/
│   ├── lbm_serial.c        # Phiên bản tuần tự 
│   └── lbm_mpi.c           # Phiên bản MPI 
├── config.txt              # Tệp cấu hình tham số mô phỏng
├── visualize_results.py    # Script vẽ trường vận tốc và mật độ
├── plot_performance.py     # Script vẽ biểu đồ so sánh hiệu năng
└── benchmark.sh            # Script tự động chạy toàn bộ quy trình test
```

### Các hàm chức năng chính
- `DocConfig()`: Đọc tham số `nx`, `ny`, `nsteps`, `omega`, `u0` từ tệp cấu hình. Trường hợp tệp thiếu hoặc lỗi, chương trình sẽ báo lỗi và dừng.
- `Collision()`: Tính toán va chạm dựa trên mô hình BGK tại các điểm lưới cục bộ.
- `Streaming()`: Lan truyền các hạt sang các điểm lưới lân cận (sử dụng mô hình Pull).
- `TraoDoiGhost()` (MPI): Thực hiện trao đổi dữ liệu biên giữa các Rank.
- `GhiKetQua()`: Xuất dữ liệu trường vận tốc và mật độ ra tệp `lbm_result.dat`.

## 3. Dữ liệu Test và các Script liên quan

### Tệp cấu hình (`config.txt`)
Mọi tham số được quản lý tập trung, không cần biên dịch lại code khi thay đổi:
```text
nx=256
ny=64
nsteps=10000
omega=1.0
u0=0.1
```

### Hướng dẫn Thực thi và Chạy Test

Người dùng có 2 tùy chọn để thực hiện mô phỏng và kiểm tra hiệu năng:

#### Tùy chọn 1: Sử dụng Script tự động (`benchmark.sh`)
Đây là cách nhanh nhất để chạy toàn bộ quy trình từ biên dịch đến vẽ biểu đồ:
```bash
chmod +x benchmark.sh
./benchmark.sh
```
*Script sẽ tự động: Biên dịch -> Chạy P=1,2,4,8 -> Vẽ biểu đồ -> Xuất thống kê.*

#### Tùy chọn 2: Thực hiện thủ công từng bước
Nếu muốn kiểm soát từng giai đoạn, sử dụng các lệnh sau:

**1. Biên dịch mã nguồn:**
```bash
gcc -O2 -o lbm_serial src/lbm_serial.c -lm
mpicc -O2 -o lbm_mpi src/lbm_mpi.c -lm
```

**2. Chạy mô phỏng:**
- **Chạy Serial:** `./lbm_serial` (mặc định dùng `config.txt`) hoặc `./lbm_serial -c config.txt`
- **Chạy MPI (ví dụ 4 tiến trình):** `mpirun -np 4 ./lbm_mpi` hoặc `mpirun -np 4 ./lbm_mpi -c config.txt`

**3. Xử lý dữ liệu và hình ảnh:**
```bash
python3 visualize_results.py  # Tạo ảnh trường vận tốc/mật độ
python3 plot_performance.py   # Cập nhật biểu đồ hiệu năng
```

### Script xử lý kết quả
- `visualize_results.py`: Tạo các ảnh `velocity_field.png`, `velocity_profile.png` (đã lược bỏ biểu đồ sai số), và `density_field.png`.
- `plot_performance.py`: Thu thập thời gian thực thi để vẽ `performance_charts.png` (Speedup, MLUPS, Efficiency).

## 4. Kết quả Thực nghiệm

### Hiệu năng (Lưới 256x64, 10,000 bước)
| Số tiến trình (P) | Thời gian (s) | MLUPS | Speedup | Efficiency |
|-------------------|---------------|-------|---------|------------|
| 1 (Serial)        | 1.594         | 102.79| 1.00×   | 100%       |
| 2                 | 0.950         | 172.55| 1.68×   | 84%        |
| 4                 | 0.691         | 236.96| 2.31×   | 58%        |
| 8                 | 1.111         | 147.48| 1.43×   | 18%        |

*Nhận xét: Hiệu năng đạt đỉnh tại 4 core thực. Tại 8 core, hiệu năng sụt giảm do tranh chấp CPU và overhead giao tiếp vượt quá khối lượng tính toán cục bộ.*

### Tính đúng đắn phục vụ so sánh (Checksum)
| Đại lượng | Serial (P=1) | MPI (P=8) | Sai lệch |
|-----------|--------------|-----------|----------|
| Sum(rho)  | 16383.9999   | 16380.9514| ~0.018%  |
| Sum(ux)   | 1527.0476    | 1527.1667 | ~0.007%  |

### Kết quả vật lý
- **Profile vận tốc**: Trùng khớp tối ưu với đường cong lý thuyết Poiseuille (RMSE = 0.0006).
- **Tính ổn định**: Trường mật độ ổn định quanh giá trị 1.0, không xảy ra hiện tượng phân kỳ tại các biên MPI.

*Báo cáo được trình bày theo yêu cầu thực nghiệm dự án HPC.*
