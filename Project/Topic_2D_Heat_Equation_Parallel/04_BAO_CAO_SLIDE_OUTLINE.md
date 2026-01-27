Khung báo cáo (10 trang)

1. Giới thiệu
- Bài toán nhiệt 2D và ứng dụng.
- Lý do chọn và mục tiêu đề tài.

2. Mô hình toán
- Phương trình, điều kiện ban đầu, điều kiện biên.
- Ý nghĩa các tham số.

3. Phương pháp số
- Lưới không gian, thời gian.
- FDM 5 điểm, công thức cập nhật.
- Điều kiện ổn định (CFL).

4. Thiết kế song song
- Chia miền 1D theo hàng.
- Hàng ghost và trao đổi biên.
- Phân tích lượng giao tiếp và đồng bộ.

5. Cài đặt
- C tuần tự, MPI.
- Cấu trúc dữ liệu, vòng lặp.
- Các tham số đầu vào.

6. Thực nghiệm
- Môi trường chạy (máy, số core).
- Tham số n, delta_t, số bước K.
- Kết quả thời gian, speedup, hiệu năng.

7. Đánh giá
- Độ chính xác (kiểm tra ổn định, tính đối xứng).
- Phân tích hiệu năng và nhận xét.

8. Kết luận
- Tổng kết và hướng mở rộng (OpenMP/CUDA, 3D, implicit).

Khung Slide thuyết trình (15 phút)

1. Tiêu đề + Thành viên thực hiện
2. Bài toán & Ý nghĩa thực tiễn
3. Mô hình toán (PDE + Điều kiện biên)
4. Phương pháp FDM 5 điểm + Điều kiện ổn định
5. Song song hóa bằng MPI: Chia miền + Ghost rows
6. Sơ đồ giao tiếp giữa các tiến trình
7. Tổng quan cài đặt mã nguồn
8. Thực nghiệm (các kịch bản tham số)
9. Kết quả speedup & hiệu năng
10. Kết luận + Hướng phát triển
