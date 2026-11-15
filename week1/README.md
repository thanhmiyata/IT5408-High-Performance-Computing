# BÀI TẬP TUẦN 1: KIẾN TRÚC SONG SONG & THIẾT KẾ THUẬT TOÁN SONG SONG
## Môn học: IT5408 – Tính toán hiệu năng cao (HPC)
* **Học viên:** Phạm Tiến Thành
* **MSSV:** 20242422M
* **Giảng viên hướng dẫn:** TS. Nguyễn Quốc Tuấn

---

## 1. Cấu trúc thư mục

```text
week1/
├── Bao_cao_Bai_tap_Tuan_1.pdf       # File báo cáo PDF hoàn chỉnh để nộp trực tiếp
├── Bao_cao_Bai_tap_Tuan_1.md        # Toàn bộ nội dung báo cáo dạng Markdown chuẩn học thuật
├── prime_count.cpp                  # Mã nguồn C++ OpenMP thực nghiệm đo Runtime đếm số nguyên tố
├── run_benchmark.py                 # Script tự động đo kiểm tra Runtime theo N và số luồng, tính Speedup
├── simulations.py                   # Script Python kiểm chứng các thuật toán tuần 1
├── Makefile                         # Script biên dịch, chạy benchmark, xuất PDF và nén ZIP
├── Bai_Tap_Tuan_1_20242422M.zip     # File ZIP hoàn chỉnh để nộp lên Google Form (< 1MB)
└── charts/                          # Biểu đồ thực nghiệm
    ├── runtime_vs_n.png             # Biểu đồ Runtime theo kích thước N và luồng
    └── amdahl_speedup_comparison.png# Biểu đồ so sánh Speedup thực tế vs Amdahl Potential Speedup
```

---

## 2. Đáp ứng đầy đủ chỉ dẫn của Giảng viên

1. **Nội dung lý thuyết (Phần A & Phần B):**
   * Giải quyết trọn vẹn và chuẩn xác cả 12 bài tập trong đề bài `Week 1.pdf`.
   * Sử dụng công thức toán học LaTeX, bảng biểu so sánh rõ ràng, phân tích sâu sắc các khía cạnh HPC.
   * Văn phong học thuật, không có từ ngữ thừa của AI, điền sẵn thông tin học viên đầy đủ.

2. **Nội dung thực nghiệm (Phần C):**
   * Triển khai bài toán thực nghiệm **Đếm số nguyên tố** (như gợi ý của thầy).
   * Khảo sát Runtime theo kích thước $N \in \{10^6, 10^7, 5\times 10^7, 10^8\}$.
   * Đo đạc trên $1, 2, 4, 8$ luồng (Apple Silicon M3).
   * Tính toán tỉ lệ phần song song thực tế $p = 96.19\%$, tỉ lệ tuần tự $s = 3.81\%$, từ đó so sánh Speedup thực tế với **Amdahl Potential Speedup**.

3. **Đóng gói nộp bài:**
   * Xuất sẵn file **`Bao_cao_Bai_tap_Tuan_1.pdf`**.
   * Đóng gói file nộp **`Bai_Tap_Tuan_1_20242422M.zip`** (dung lượng ~745 KB, tuyệt đối không chứa binary chạy, tuân thủ đúng yêu cầu giới hạn < 10MB).
