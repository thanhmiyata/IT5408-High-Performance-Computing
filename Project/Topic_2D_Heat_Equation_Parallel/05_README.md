README - 2D Heat Equation (MPI)

1. Mô tả
- Giải phương trình nhiệt 2D bằng phương pháp FDM explicit.
- Song song hóa bằng thư viện MPI (chia miền 1D theo hàng).

2. Cấu trúc thư mục xu đề xuất
- src/diffusion_serial.c: Mã nguồn tuần tự.
- src/diffusion_mpi.c: Mã nguồn song song MPI.
- data/: Thư mục chứa các file cấu hình đầu vào.
- report.pdf: Báo cáo chi tiết.
- slides.pdf: Slide thuyết trình.

3. Dữ liệu đầu vào (Input data)
- File cấu hình sử dụng định dạng key=value, ví dụ:
  n=100
  K=500
  c=0.1
  left=10
  right=40
  bottom=30
  top=50
  init=0
- Các bộ dữ liệu có sẵn:
  data/case1.txt
  data/case2.txt

4. Biên dịch (Compilation)
- Tuần tự (Serial):
  gcc -O2 -o diffusion_serial src/diffusion_serial.c -lm
- Song song (MPI):
  mpicc -O2 -o diffusion_mpi src/diffusion_mpi.c -lm

5. Chạy chương trình (Execution)
- Tuần tự:
  ./diffusion_serial data/case1.txt
- Song song:
  mpirun -np 4 ./diffusion_mpi data/case1.txt

Giải thích tham số:
- n: số điểm nội bộ mỗi chiều (tổng lưới là n+2).
- K: số bước thời gian.
- c: hệ số khuếch tán.
- left/right/bottom/top: nhiệt độ tại các biên tương ứng.
- init: giá trị nhiệt độ ban đầu bên trong lưới.

6. Kiểm tra kết quả
- Kiểm tra điều kiện ổn định: delta_t <= delta_s^2 / (2c).
- In ra giá trị u tại một vài điểm đại diện để so sánh.
- So sánh giá trị tổng (checksum) giữa bản tuần tự và song song.

7. Đo lường hiệu năng
- Lấy thời gian chạy bằng hàm MPI_Wtime (tại rank 0).
- Speedup = T1 / Tp, Efficiency = Speedup / p.
