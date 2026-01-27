Thiết kế song song MPI (chia miền 1D)

1. Ý tưởng chia miền
- Chia lưới theo chiều x (chia theo hàng).
- Mỗi tiến trình (process) xử lý một khối hàng liên tiếp.
- Cần 2 hàng ảo (ghost rows) để cập nhật stencil 5 điểm.

2. Kích thước cục bộ
Giả sử số hàng nội bộ là n.
Số tiến trình = P.
Chia đó:
rows_per_proc = n / P, xử lý phần dư nếu n % P != 0.

Mỗi tiến trình có:
- local_n = số hàng nội bộ mà tiến trình quản lý.
- Mảng cục bộ: (local_n + 2) x (n + 2)
  (bao gồm 2 hàng ghost).

3. Trao đổi biên (ghost rows)
Trước mỗi bước thời gian:
- Nếu rank > 0: gửi hàng đầu nội bộ (1) cho rank-1,
  nhận ghost trên (0) từ rank-1.
- Nếu rank < P-1: gửi hàng cuối nội bộ (local_n) cho rank+1,
  nhận ghost dưới (local_n+1) từ rank+1.

Có thể dùng MPI_Sendrecv hoặc MPI_Isend/Irecv.

4. Vòng lặp cập nhật
for k in 1..K:
  trao đổi các hàng ghost
  for i in 1..local_n:
    for j in 1..n:
      cập nhật u_new[i][j] theo công thức FDM
  swap(u_old, u_new)

5. Điều kiện biên
- Hàng biên toàn cục: i = 0 và i = n+1.
- Chỉ tiến trình đầu/cuối mới gán nhiệt độ biên tương ứng.
- Các cột biên (j=0, j=n+1) gán tương ứng trên tất cả tiến trình.

6. Gom kết quả
Sau khi xong, gom các khối về rank 0:
Sử dụng MPI_Gatherv trên các hàng nội bộ.

7. Đo hiệu năng
- Đo thời gian bằng MPI_Wtime (rank 0).
- Tính speedup và hiệu năng (efficiency).

8. Lời khuyên
- Nếu n nhỏ, chi phí liên lạc lớn, speedup kém.
- Nếu n lớn, speedup sẽ tốt hơn.
