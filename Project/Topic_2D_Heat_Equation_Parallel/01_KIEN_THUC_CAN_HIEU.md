Nội dung cần hiểu (giải thích chi tiết):

1. Phương trình nhiệt 2D
- Đối tượng: nhiệt độ u(x, y, t) trên miền vuông [0,1]x[0,1] theo thời gian.
- Phương trình: u_t = c * (u_xx + u_yy), c > 0 là hệ số khuếch tán.
- Ý nghĩa: nhiệt lan truyền từ nơi nhiệt cao sang thấp theo thời gian.

2. Điều kiện biên và điều kiện ban đầu
- Điều kiện ban đầu: u(0, x, y) = f(x, y). Cho biết nhiệt độ ban đầu.
- Điều kiện biên: u(t, 0, y) = alpha0(y), u(t, 1, y) = alpha1(y),
  u(t, x, 0) = beta0(x), u(t, x, 1) = beta1(x).
- Trong lập trình, biên được gán cố định và không cập nhật theo thời gian.

3. Rời rạc hóa miền (lưới không gian và thời gian)
- Chia đoạn [0,1] thành n+1 khoảng, bước lưới delta_s = 1/(n+1).
- Các điểm lưới: x_i = i * delta_s, y_j = j * delta_s, i,j = 0..n+1.
- Thời gian: t_k = k * delta_t, k = 0,1,2,...

4. Sai phân hữu hạn (Finite Difference Method - FDM)
- Đạo hàm theo thời gian:
  u_t ~ (u^{k+1}_{i,j} - u^k_{i,j}) / delta_t.
- Đạo hàm bậc hai theo không gian:
  u_xx ~ (u^k_{i+1,j} - 2u^k_{i,j} + u^k_{i-1,j}) / delta_s^2
  u_yy ~ (u^k_{i,j+1} - 2u^k_{i,j} + u^k_{i,j-1}) / delta_s^2
- Thay vào phương trình nhiệt ta có công thức cập nhật (explicit):
  u^{k+1}_{i,j} = u^k_{i,j}
    + r * (u^k_{i+1,j} + u^k_{i-1,j} + u^k_{i,j+1} + u^k_{i,j-1} - 4u^k_{i,j})
  với r = c * delta_t / delta_s^2.

5. Điều kiện ổn định (CFL)
- FDM explicit ổn định nếu delta_t <= (delta_s^2) / (2c).
- Nếu vi phạm, kết quả sẽ dao động và phân kỳ (diverge).
- Trong mã nguồn, chọn delta_t theo công thức trên để an toàn.

6. Kiểm chứng kết quả
- Kiểm tra tính đối xứng: nếu điều kiện biên đối xứng thì nghiệm đối xứng.
- Kiểm tra đạo hàm: u không bị lỗi số (NaN).
- So sánh với nghiệm tiêu chuẩn: bài toán với điều kiện biên đơn giản có thể có nghiệm ổn định (steady-state).
- Kiểm tra hội tụ: khi tăng n (lưới tinh hơn), kết quả ổn định hơn.

7. Hiểu thuật toán tuần tự
- Tuần tự: vòng lặp theo thời gian, mỗi bước tính toán toàn bộ lưới.
- Mỗi điểm (i,j) phụ thuộc 4 hàng xóm ở thời điểm k.
- Tính xong u_new cho toàn bộ, rồi hoán đổi u_old/u_new.

8. Lý do có thể song song hóa
- Cập nhật u^{k+1}_{i,j} chỉ cần 4 điểm láng giềng (i+1,i-1,j+1,j-1).
- Do đó có thể chia miền theo hàng (y) hay cột (x).
- Chỉ cần trao đổi biên (ghost rows/cols) với tiến trình kế bên.

9. MPI cần biết
- MPI_Init, MPI_Comm_rank, MPI_Comm_size.
- MPI_Sendrecv hoặc MPI_Isend/Irecv để trao đổi biên.
- MPI_Barrier để đồng bộ.
- Đo thời gian bằng MPI_Wtime.

10. Đo hiệu năng
- Time1: thời gian chạy với 1 tiến trình.
- TimeP: thời gian chạy với P tiến trình.
- Speedup = Time1 / TimeP.
- Efficiency = Speedup / P.

11. Định dạng dữ liệu
- Mảng 2D (n+2 x n+2) gồm biên + nội bộ.
- u_old và u_new, hoán đổi mỗi bước thời gian.
- Chỉ cập nhật điểm nội bộ (1..n).
