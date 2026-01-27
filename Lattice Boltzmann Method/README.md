# Lattice Boltzmann Method (LBM) - Mô phỏng Dòng chảy 2D

## Giới thiệu Đề tài

Đề tài này nghiên cứu **phương pháp Lattice Boltzmann (LBM)** để mô phỏng dòng chảy chất lỏng trong không gian 2 chiều, đặc biệt tập trung vào việc **song song hóa** thuật toán trên các hệ thống tính toán hiệu năng cao.

### Tại sao chọn LBM?

1. **Tính song song tự nhiên**: Mỗi điểm lưới chỉ tương tác với các điểm lân cận → dễ dàng phân chia miền
2. **Ứng dụng thực tế**: Khí động học ô tô, thiết kế cánh máy bay, mô phỏng dòng chảy trong vi mạch
3. **Hiệu quả trên GPU**: Cấu trúc dữ liệu đều đặn, phù hợp với kiến trúc SIMD
4. **Xử lý biên phức tạp**: Dễ dàng mô phỏng các vật cản hình dạng bất kỳ

## Bài toán Cụ thể

**Mô phỏng dòng chảy Poiseuille (Dòng chảy giữa 2 tấm phẳng song song)**

- Miền tính toán: Lx × Ly (ví dụ: 512 × 128 điểm lưới)
- Điều kiện biên:
  - Trên/Dưới: Tường cứng (bounce-back)
  - Trái: Vận tốc đầu vào u₀
  - Phải: Áp suất đầu ra
- Số Reynolds: Re = 100-1000

## Cấu trúc Thư mục

```
Lattice Boltzmann Method/
├── README.md              # File này
├── THEORY.md              # Lý thuyết chi tiết về LBM
├── src/
│   ├── lbm_sequential.c   # Phiên bản tuần tự (baseline)
│   ├── lbm_mpi.c          # Phiên bản MPI
│   ├── lbm_openmp.c       # Phiên bản OpenMP
│   └── lbm_hybrid.c       # Phiên bản MPI + OpenMP
├── include/
│   └── lbm.h              # Header chung
├── scripts/
│   ├── visualize.py       # Script Python để vẽ kết quả
│   └── benchmark.sh       # Script đo hiệu năng
├── Makefile
└── results/               # Thư mục chứa kết quả
```

## Yêu cầu Hệ thống

- **Compiler**: GCC 9.0+ hoặc Intel ICC
- **MPI**: OpenMPI 4.0+ hoặc MPICH 3.3+
- **OpenMP**: Hỗ trợ OpenMP 4.5+
- **Python** (tùy chọn): Python 3.8+ với matplotlib, numpy để visualize

## Biên dịch

```bash
# Phiên bản tuần tự
make sequential

# Phiên bản MPI
make mpi

# Phiên bản OpenMP
make openmp

# Phiên bản Hybrid
make hybrid

# Biên dịch tất cả
make all
```

## Chạy thử nghiệm

```bash
# Tuần tự
./bin/lbm_sequential 512 128 10000

# MPI với 4 tiến trình
mpirun -np 4 ./bin/lbm_mpi 512 128 10000

# OpenMP với 8 luồng
export OMP_NUM_THREADS=8
./bin/lbm_openmp 512 128 10000

# Hybrid: 2 tiến trình MPI × 4 luồng OpenMP
export OMP_NUM_THREADS=4
mpirun -np 2 ./bin/lbm_hybrid 512 128 10000
```

## Phân tích Hiệu năng

### Các chỉ số đo lường:

1. **MLUPS** (Million Lattice Updates Per Second): Số triệu điểm lưới được cập nhật mỗi giây
2. **Speedup**: T_sequential / T_parallel
3. **Efficiency**: Speedup / Số_processors
4. **Băng thông bộ nhớ**: GB/s

### Kết quả mong đợi:

| Phiên bản | Lưới 512×128 | Lưới 2048×512 | Speedup (8 cores) |
|-----------|--------------|---------------|-------------------|
| Sequential| ~50 MLUPS    | ~45 MLUPS     | 1.0×              |
| OpenMP    | ~350 MLUPS   | ~320 MLUPS    | 7.1×              |
| MPI       | ~380 MLUPS   | ~340 MLUPS    | 7.6×              |
| Hybrid    | ~400 MLUPS   | ~360 MLUPS    | 8.0×              |

## Tài liệu Tham khảo

1. **Sách**: "The Lattice Boltzmann Method: Principles and Practice" - Krüger et al. (2017)
2. **Paper**: "Lattice Boltzmann Method for Fluid Flows" - Chen & Doolen (1998)
3. **Code mẫu**: [palabos.org](http://www.palabos.org/)
4. **Tutorial**: [lbmethod.org](http://www.lbmethod.org/)

## Tác giả

- **Sinh viên**: Krizpham
- **Môn học**: IT5408 - Tính toán Hiệu năng Cao
- **Trường**: [Tên trường của bạn]
- **Năm**: 2026

## License

MIT License - Tự do sử dụng cho mục đích học tập và nghiên cứu.
