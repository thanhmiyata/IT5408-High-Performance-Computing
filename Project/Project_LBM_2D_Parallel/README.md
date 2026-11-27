# 🚀 Mô phỏng Dòng chảy 2D bằng Lattice Boltzmann Method (LBM D2Q9) và MPI

**Học phần:** IT5408 - Tính toán hiệu năng cao  
**Giảng viên hướng dẫn:** TS. Nguyễn Quốc Tuấn  
**Đề tài:** 2D Hydrodynamic Flow Simulation - Parallel Numerical Solution (MPI)  
**GitLab Repository:** `https://gitlab.nqtuan0192.me/<username>/HPC_20261B_Group<Number>_LBM_MPI`

---

## 👥 Danh sách thành viên nhóm

| STT | Họ và tên | Mã học viên / MSSV | Vai trò & Đóng góp |
|:---:|:---|:---:|:---|
| 1 | **Phạm Tiến Thành** | **20242422M** | Trưởng nhóm, Cài đặt song song MPI, Benchmark hiệu năng |
| 2 | **Trần Thị Hạnh** | **20242058M** | Phân tích mô hình toán D2Q9, Điều kiện biên vật lý |
| 3 | **Bùi Minh Đức** | **20242056M** | Cài đặt thuật toán tuần tự C99, Kiểm thử Checksum |
| 4 | **Nguyễn Thị Huê** | **20252584M** | Xử lý dữ liệu, Script visualize trường vận tốc/mật độ |
| 5 | **Nguyễn Đình Hiếu** | **20261223M** | Thiết kế slide thuyết trình Beamer, Phân tích Scaling |
| 6 | **Đỗ Minh Hiếu** | **20242031M** | Biên soạn báo cáo kỹ thuật, Phân tích Overhead giao tiếp |

---

## 📋 Tổng quan Đề tài

Dự án cài đặt phương pháp **Lattice Boltzmann Method (LBM)** để mô phỏng dòng chảy chất lưu 2 chiều trong kênh dẫn phẳng (Poiseuille Flow), so sánh giữa phiên bản tuần tự (Serial C99) và phiên bản song song hóa phân tán (**MPI**):
- **Mô hình lưới vi mô:** D2Q9 (2D, 9 vector vận tốc rời rạc).
- **Mô hình va chạm:** LBGK (Bhatnagar-Gross-Krook) với Single Relaxation Time $\tau$.
- **Chiến lược song song:** Chia miền 1 chiều (**1D Domain Decomposition**) theo các cột trục $X$.
- **Giao tiếp liên tiến trình:** Trao đổi các cột biên ảo (**Ghost Columns**) sử dụng kỹ thuật không chặn (**Non-blocking**: `MPI_Irecv` và `MPI_Isend` kết hợp `MPI_Waitall`).
- **Thu thập dữ liệu:** Sử dụng `MPI_Gatherv` để gom toàn bộ dữ liệu trường vận tốc/mật độ về Rank 0 ghi file, và `MPI_Reduce` đối soát Checksum.
- **Tính độc lập:** Tự thiết kế và cài đặt từng dòng code C (không sử dụng bất kỳ thư viện ngoài nào).

---

## 📁 Cấu trúc thư mục dự án

```text
Project_LBM_2D_Parallel/
├── src/
│   ├── lbm_serial.c           # Mã nguồn C tuần tự (Baseline)
│   └── lbm_mpi.c              # Mã nguồn C song song hóa bằng MPI
├── data/
│   ├── config_small.txt       # Cấu hình lưới nhỏ (128×32)
│   ├── config_medium.txt      # Cấu hình mặc định (256×64)
│   └── config_large.txt       # Cấu hình lưới lớn (512×128)
├── Makefile                   # Tự động hóa biên dịch và thực thi
├── benchmark.sh               # Script tự động chạy benchmark toàn diện (P=1,2,4,8)
├── plot_performance.py        # Vẽ biểu đồ Speedup, MLUPS, Efficiency
├── visualize_results.py       # Trực quan hóa trường vận tốc, vector, mật độ
├── config.txt                 # File cấu hình tham số hiện hành
├── Slides 0402.pdf            # Slide thuyết trình Beamer (29 trang)
├── Bao_cao_HPC_LBM.md         # Báo cáo kỹ thuật chi tiết (~10 trang)
├── Bao_cao_HPC_LBM.pdf        # Báo cáo kỹ thuật định dạng PDF chuẩn nộp bài
└── *.png                      # Các biểu đồ hiệu năng và ảnh mô phỏng chất lưu
```

---

## 🛠️ Yêu cầu môi trường & Cài đặt

1. **Trình biên dịch C:** GCC (hỗ trợ C99 trở lên) hoặc Clang.
2. **MPI Runtime:** OpenMPI hoặc MPICH (`mpicc`, `mpirun`).
3. **Python 3:** Thư viện `numpy`, `matplotlib` (để visualize kết quả).

---

## 🚀 Hướng dẫn Biên dịch và Chạy

### Cách 1: Sử dụng `Makefile` (Khuyên dùng)

```bash
# 1. Biên dịch cả 2 phiên bản Serial và MPI
make

# 2. Chạy bản tuần tự
make run_serial

# 3. Chạy bản song song MPI với 4 tiến trình (mặc định)
make run_mpi NP=4

# 4. Chạy demo kèm UI trực quan xem luồng nước chảy theo thời gian (lưu mỗi 200 bước)
make demo NP=4

# 5. Trực quan hóa kết quả mô phỏng và cập nhật biểu đồ
make plot

# 6. Dọn dẹp file thực thi và file dữ liệu tạm
make clean
```

### Cách 2: Sử dụng lệnh trực tiếp

```bash
# Biên dịch
gcc -O2 -Wall -o lbm_serial src/lbm_serial.c -lm
mpicc -O2 -Wall -o lbm_mpi src/lbm_mpi.c -lm

# Chạy bản tuần tự
./lbm_serial -c config.txt

# Chạy MPI với P tiến trình (ví dụ 4 tiến trình)
mpirun -np 4 ./lbm_mpi -c config.txt

# Chạy kèm tham số -f / --frames để xuất chuỗi frame cho UI xem dòng chảy
mpirun -np 4 ./lbm_mpi -c config.txt -f 200

# Trực quan hóa dữ liệu (ảnh tĩnh + dữ liệu UI web + GIF nếu thêm --gif)
python3 visualize_results.py --demo
python3 plot_performance.py
```

### Cách 3: Chạy toàn bộ quy trình kiểm thử tự động
```bash
chmod +x benchmark.sh
./benchmark.sh
```

---

## 📊 Kết quả Thực nghiệm & Hiệu năng

### Thử nghiệm trên lưới $256 \times 64$, 10,000 bước thời gian ($\omega = 1.0, u_0 = 0.1$):

| Số tiến trình ($P$) | Thời gian thực thi ($s$) | MLUPS | Speedup ($S_p$) | Efficiency ($E_p$) |
|:---:|:---:|:---:|:---:|:---:|
| **1 (Serial)** | **1.584** | **103.45** | **1.00×** | **100.0%** |
| **2** | **0.906** | **180.87** | **1.75×** | **87.4%** |
| **4** | **0.507** | **323.26** | **3.12×** | **78.1%** |
| **8** | **1.355** | **120.92** | **1.17×** | **14.6%** |

* **MLUPS (Million Lattice Updates Per Second):** Số triệu điểm lưới được cập nhật mỗi giây, đại lượng đo hiệu năng tiêu chuẩn của mô phỏng LBM:
  $$\text{MLUPS} = \frac{N_x \times N_y \times N_{steps}}{T \times 10^6}$$
* **Nhận xét hiệu năng:**
  * Hiệu năng đạt đỉnh tại $P = 4$ cores thực tế, tăng tốc độ mô phỏng gấp **3.12 lần** so với bản tuần tự, đưa MLUPS từ 103.45 lên **323.26** (hiệu suất 78.1%).
  * Khi tăng lên $P = 8$, hiệu năng giảm do overhead truyền thông biên giữa các rank bắt đầu vượt qua khối lượng tính toán cục bộ (mỗi tiến trình chỉ còn $256/8 = 32$ cột) kết hợp hiện tượng tranh chấp tài nguyên trên chip hybrid.

---

## 🔍 Kiểm chứng Tính đúng đắn (Checksum & Vật lý)

### 1. Đối chiếu Checksum giữa Serial và MPI ($P=4$):
* $\sum \rho$ (Serial): `16383.9999999909` $\longleftrightarrow$ $\sum \rho$ (MPI): `16383.9999999910` (Khớp tuyệt đối, sai lệch $< 10^{-12}\%$)
* $\sum u_x$ (Serial): `1527.0476269322` $\longleftrightarrow$ $\sum u_x$ (MPI): `1527.0476269322` (Khớp tuyệt đối 10 chữ số thập phân)
* Khối lượng chất lưu được bảo toàn hoàn hảo qua 10,000 bước lặp và xuyên suốt các biên chia miền MPI.

### 2. Kết quả vật lý:
* Profile vận tốc dọc theo trục $Y$ khớp chính xác với đường cong Parabol giải tích của dòng Poiseuille:
  $$u_x(y) = 4 u_{max} \frac{y(H - y)}{H^2}$$
* Hình ảnh trường vận tốc (`velocity_field.png`) và profile vận tốc (`velocity_profile.png`) thể hiện sự ổn định và liên tục qua các biên chia miền MPI.
