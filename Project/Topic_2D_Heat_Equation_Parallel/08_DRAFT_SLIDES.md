# SLIDE DỰ THẢO: GIẢI PHƯƠNG TRÌNH NHIỆT 2D (MPI)
**Thời lượng dự kiến:** 15 phút

---

## Slide 1: Tiêu đề
- **Tên đề tài:** 2D Heat Equation - Parallel Numerical Solution with MPI.
- **Môn học:** IT5408 - Tính toán hiệu năng cao.
- **Thành viên thực hiện:** [Tên của bạn].

---

## Slide 2: Đặt vấn đề & Ứng dụng
- Bài toán: Mô phỏng phân bố nhiệt độ trên vật thể 2D.
- Ứng dụng:
    - Làm mát linh kiện điện tử.
    - Dự báo lan truyền nhiệt trong vật liệu.
    - Cơ sở cho các bài toán động lực học chất lưu (CFD).

---

## Slide 3: Mô hình toán học
- Phương trình: $u_t = c \cdot (u_{xx} + u_{yy})$.
- Miền tính toán: Lưới vuông đơn vị.
- Điều kiện biên Dirichlet: Nhiệt độ 4 cạnh cố định.
- Điều kiện ban đầu: Phân bố nhiệt tại $t=0$.

---

## Slide 4: Phương pháp số (FDM)
- Rời rạc hóa không gian và thời gian.
- **Stencil 5 điểm:** Một điểm phụ thuộc vào chính nó và 4 điểm lân cận (trên, dưới, trái, phải).
- Điều kiện ổn định (CFL): $\Delta t$ phải đủ nhỏ để tránh phân kỳ kết quả.

---

## Slide 5: Chiến lược song song hóa
- Sử dụng mô hình truyền tin MPI.
- **Chia miền 1D:** Chia theo khối hàng (Row-wise decomposition).
- Ưu điểm: Đơn giản trong việc lập trình và quản lý bộ nhớ.

---

## Slide 6: Vùng đệm ảo (Ghost Rows)
- Thách thức: Điểm biên của khối này cần dữ liệu từ khối kia.
- Giải pháp: Thêm 2 hàng ảo (Ghost Rows) ở phía trên và dưới mỗi khối dữ liệu cục bộ.
- Hình ảnh minh họa: [Sơ đồ trao đổi biên giữa Rank $i$ và Rank $i \pm 1$].

---

## Slide 7: Cơ chế giao tiếp MPI
- Sử dụng `MPI_Sendrecv`: Tránh tình trạng Deadlock (khóa chết) khi các tiến trình gửi tin cho nhau đồng thời.
- Quy trình:
    1. Rank $i$ gửi hàng $1$ lên biên trên của Rank $i-1$.
    2. Rank $i$ gửi hàng $local\_n$ xuống biên dưới của Rank $i+1$.

---

## Slide 8: Cài đặt và Luồng chương trình
- Ngôn ngữ: C + MPI.
- Các bước chính:
    1. Đọc file config.
    2. Phân chia dữ liệu (MPI_Scatter/Gatherv).
    3. Vòng lặp tính toán & Hoán đổi biên.
    4. Tổng hợp và xuất kết quả.

---

## Slide 9: Thực nghiệm (Kết quả & Biểu đồ)
- Thông số lưới: $1000 \times 1000$ (hoặc lớn hơn).
- Biểu đồ thời gian chạy theo số Core.
- Biểu đồ Speedup và Efficiency.
- *(Để trống số liệu để điền sau)*.

---

## Slide 10: Phân tích hiệu năng
- Tại sao Speedup không đạt mức lý tưởng ($P$ lần)?
    - Chi phí khởi tạo MPI.
    - Thời gian trễ đường truyền (Networking Latency) khi trao đổi Ghost Rows.
    - Phần code tuần tự (Đọc file, ghi kết quả) theo định luật Amdahl.

---

## Slide 11: Kết luận & Hướng phát triển
- Đạt được: Giải bài toán chính xác, hiểu sâu về phân chia dữ liệu và MPI.
- Dự kiến mở rộng:
    - Chia miền 2D (Checkerboard) để giảm lượng dữ liệu trao đổi.
    - Kết hợp MPI + OpenMP (Hybrid) để tận dụng đa nhân trên mỗi Node.

---

## Slide 12: Q&A
- Cảm ơn thầy và các bạn đã lắng nghe.
- Câu hỏi thảo luận.
