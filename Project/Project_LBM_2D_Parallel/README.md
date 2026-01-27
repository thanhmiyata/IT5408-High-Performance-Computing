# 🚀 Lattice Boltzmann Method - Mô phỏng Dòng chảy 2D Song song

## 📋 Tổng quan

Dự án này cài đặt **Lattice Boltzmann Method (LBM)** để mô phỏng dòng chảy chất lưu 2D, với phiên bản tuần tự và song song MPI.

- **Mô hình:** D2Q9 (2D, 9 hướng vận tốc)
- **Bài toán:** Dòng chảy Poiseuille
- **Song song hóa:** MPI với chia miền 1D

## 📁 Cấu trúc dự án

```
Project_LBM_2D_Parallel/
├── src/
│   ├── lbm_serial.c           # Code tuần tự
│   └── lbm_mpi.c              # Code song song MPI
├── data/
│   ├── config_small.txt       # Cấu hình test nhỏ (128×32)
│   └── config_large.txt       # Cấu hình test lớn (512×128)
├── benchmark.sh               # Script chạy benchmark tự động
├── plot_performance.py        # Vẽ biểu đồ hiệu năng
├── visualize_results.py       # Visualize kết quả mô phỏng
├── 00-08_*.md                 # Tài liệu chi tiết
└── *.png                      # Biểu đồ kết quả
```

## 🚀 Hướng dẫn sử dụng nhanh

### 1. Biên dịch
```bash
# Serial
gcc -O2 -o lbm_serial src/lbm_serial.c -lm

# MPI
mpicc -O2 -o lbm_mpi src/lbm_mpi.c -lm
```

### 2. Chạy
```bash
# Serial
./lbm_serial

# MPI với 4 tiến trình
mpirun -np 4 ./lbm_mpi
```

### 3. Benchmark tự động
```bash
chmod +x benchmark.sh
./benchmark.sh
```

### 4. Vẽ biểu đồ
```bash
# Biểu đồ hiệu năng (Speedup, MLUPS, Efficiency)
python3 plot_performance.py

# Visualize kết quả (trường vận tốc, profile)
python3 visualize_results.py
```

## 📊 Kết quả hiệu năng

### Lưới 256×64, 10,000 bước

| P | Thời gian (s) | MLUPS | Speedup | Efficiency |
|---|---------------|-------|---------|------------|
| 1 | 1.802 | 90.94 | 1.00× | 100% |
| 2 | 1.365 | 120.00 | 1.32× | 66% |
| 4 | 0.738 | 222.13 | 2.44× | 61% |
| 8 | 0.450 | 364.09 | 4.00× | 50% |

**MLUPS:** Million Lattice Updates Per Second

## 📈 Biểu đồ

Sau khi chạy `plot_performance.py` và `visualize_results.py`, bạn sẽ có:

1. **performance_charts.png** - Speedup, MLUPS, Efficiency
2. **velocity_field.png** - Trường vận tốc (contour + vector)
3. **velocity_profile.png** - Profile vận tốc theo chiều Y
4. **density_field.png** - Trường mật độ

## 📚 Tài liệu

- `00_OVERVIEW.md` - Tổng quan đề tài
- `01_KIEN_THUC_CAN_HIEU.md` - Kiến thức về LBM
- `02_MO_HINH_TOAN.md` - Mô hình toán học D2Q9
- `03_THIET_KE_SONG_SONG_MPI.md` - Thiết kế song song
- `04_BAO_CAO_SLIDE_OUTLINE.md` - Khung báo cáo và slide
- `05_README.md` - Hướng dẫn chi tiết
- `06_TEST_DANH_GIA.md` - Hướng dẫn test
- `07_DRAFT_REPORT.md` - Báo cáo đầy đủ
- `08_DRAFT_SLIDES.md` - Slide thuyết trình

## 🔧 Tùy chỉnh tham số

Chỉnh sửa các `#define` trong file `.c`:

```c
#define NX  256      // Kích thước X
#define NY  64       // Kích thước Y
#define NSTEPS 10000 // Số bước thời gian
#define OMEGA  1.0   // Tham số thư giãn
#define U0     0.1   // Vận tốc đầu vào
```

Hoặc sử dụng file config (tính năng đang phát triển).

## 🎯 Mục tiêu đạt được

✅ Cài đặt LBM D2Q9 tuần tự và MPI  
✅ Speedup 2.44× với 4 tiến trình  
✅ MLUPS tăng từ 90.94 → 222.13  
✅ Visualize kết quả mô phỏng  
✅ Tài liệu đầy đủ (báo cáo + slide)  

## 🚧 Hướng phát triển

- [ ] Chia miền 2D để giảm giao tiếp
- [ ] Hybrid MPI + OpenMP
- [ ] GPU (CUDA/OpenCL)
- [ ] Mở rộng sang 3D (D3Q19)
- [ ] Điều kiện biên Zou-He
- [ ] Đọc tham số từ file config

## 📖 Tài liệu tham khảo

1. T. Krüger et al., "The Lattice Boltzmann Method" (2017)
2. S. Chen & G.D. Doolen, "Lattice Boltzmann Method for Fluid Flows" (1998)
3. http://www.lbmethod.org/
4. http://www.palabos.org/

## 👤 Tác giả

- **Sinh viên:** Krizpham
- **Môn học:** IT5408 - Tính toán Hiệu năng Cao
- **Năm:** 2026

## 📝 License

MIT License - Tự do sử dụng cho mục đích học tập và nghiên cứu.
