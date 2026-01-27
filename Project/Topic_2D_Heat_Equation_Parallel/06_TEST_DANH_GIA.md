Hướng dẫn chạy thử nghiệm và đánh giá

1. Mục tiêu thử nghiệm
- So sánh tính đúng đắn giữa bản tuần tự (serial) và bản song song (MPI).
- Đánh giá hiệu suất: thời gian chạy, speedup, hiệu năng (efficiency).

2. Kiểm tra tính đúng đắn
- Chạy bản tuần tự và song song với cùng một bộ tham số đầu vào.
- So sánh các giá trị:
  + u_mid (giá trị nhiệt độ tại tâm lưới).
  + sum (tổng tất cả các giá trị u trên lưới).
- Hai giá trị này cần khớp nhau hoặc có sai số cực nhỏ (do làm tròn số thực).

3. Đánh giá thời gian
Chạy thử nghiệm nhiều lần và lấy giá trị trung bình để đảm bảo khách quan:
- Bản tuần tự:
  ./diffusion_serial data/case2.txt
- Bản MPI:
  mpirun -np 2 ./diffusion_mpi data/case2.txt
  mpirun -np 4 ./diffusion_mpi data/case2.txt
  mpirun -np 8 ./diffusion_mpi data/case2.txt

Ghi chép lại thời gian (Time) từ kết quả đầu ra.

4. Tính toán Speedup và Hiệu năng
Speedup(p) = T1 / Tp
Efficiency(p) = Speedup(p) / p

5. Kỳ vọng kết quả
- Khi kích thước lưới (n) càng lớn, hiệu quả song song hóa (speedup) càng rõ rệt.
- Nếu lưới n quá nhỏ, chi phí trao đổi dữ liệu qua mạng sẽ áp đảo thời gian tính toán, dẫn đến hiệu năng giảm.