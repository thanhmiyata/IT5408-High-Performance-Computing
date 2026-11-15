# BÀI TẬP TUẦN 2: LẬP TRÌNH SONG SONG VỚI OPENMP & MPI
## Môn học: IT5408 – Tính toán hiệu năng cao (HPC)
* **Học viên:** Phạm Tiến Thành
* **MSSV:** 20242422M
* **Giảng viên hướng dẫn:** TS. Nguyễn Quốc Tuấn

---

## 1. Cấu trúc thư mục

```text
Week 2/
├── Bao_cao_Bai_tap_Tuan_2.pdf       # File báo cáo PDF hoàn chỉnh để nộp trực tiếp
├── Bao_cao_Bai_tap_Tuan_2.md        # Toàn bộ nội dung báo cáo chuẩn học thuật
├── tsp_solver.cpp                   # Bài 1: Giải TSP bằng Nhánh & Cận song song OpenMP
├── prime_count.cpp                  # Bài 2: Đếm số nguyên tố song song OpenMP
├── count_paths.cpp                  # Bài 3: Đếm số đường đi độ dài k bằng lũy thừa ma trận kề
├── matmul_serial.cpp                # Bài bổ sung: Nhân ma trận tuần tự
├── matmul_thread.cpp                # Bài bổ sung: Nhân ma trận với C++ std::thread
├── matmul_openmp.cpp                # Bài bổ sung: Nhân ma trận với OpenMP
├── matmul_mpi.cpp                   # Bài bổ sung: Nhân ma trận với OpenMPI
├── common.hpp                       # Header chung cho nhân ma trận
├── benchmark_week2.py               # Script tự động đo kiểm toàn bộ 4 bài toán và vẽ biểu đồ
├── export_pdf.py                    # Script xuất báo cáo PDF và nén ZIP
├── Makefile                         # Script build tự động
├── Bai_Tap_Tuan_2_20242422M.zip     # File ZIP nộp bài (< 1MB, không chứa binary)
├── charts/                          # Biểu đồ thực nghiệm
│   ├── prime_runtime_amdahl.png
│   ├── tsp_scaling.png
│   └── matmul_comparison.png
├── tsp*.inp                         # Dữ liệu kiểm thử TSP
└── graph*.txt                       # Dữ liệu kiểm thử đồ thị
```

---

## 2. Hướng dẫn biên dịch và chạy

```bash
cd "Week 2"

# 1. Biên dịch toàn bộ
make all

# 2. Chạy thử nghiệm nhanh các bài toán
make test

# 3. Chạy từng bài riêng biệt:
./tsp_solver tsp15.inp 4
./prime_count 10000000 4
./count_paths graph05.txt 5 4
./matmul_openmp 1024 1024 1024 4

# 4. Chạy lại benchmark toàn diện và vẽ đồ thị:
source .venv/bin/activate
python3 benchmark_week2.py

# 5. Xuất lại file PDF và đóng gói ZIP nộp bài:
python3 export_pdf.py
```
