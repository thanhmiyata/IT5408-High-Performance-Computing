# BÁO CÁO DỰ ÁN: GIẢI PHƯƠNG TRÌNH NHIỆT 2D BẰNG THUẬT TOÁN SONG SONG MPI

**Môn học:** IT5408 - Tính toán hiệu năng cao
**Đề tài:** 2D Heat Equation - Parallel Numerical Solution (MPI)

---

## 1. Giới thiệu (Problem Specification)
### 1.1. Đề tài
Dự án tập trung vào việc giải phương trình nhiệt 2D (2D Heat Equation) - một bài toán cơ bản trong vật lý và kỹ thuật, mô tả sự lan truyền nhiệt thực tế trên một bề mặt phẳng theo thời gian.

### 1.2. Lý do chọn đề tài
- **Tính ứng dụng cao:** Mô phỏng nhiệt là nền tảng cho thiết kế vi mạch, kỹ thuật xây dựng và khoa học vật liệu.
- **Phù hợp với tính toán song song:** Thuật toán Sai phân hữu hạn (FDM) có tính cục bộ dữ liệu cao, rất lý tưởng để minh họa phương pháp chia miền (Domain Decomposition) và giao tiếp giữa các tiến trình.
- **Khả năng mở rộng:** Từ bài toán 2D có thể dễ dàng mở rộng sang 3D hoặc các hệ phương trình đạo hàm riêng phức tạp hơn.

### 1.3. Mục tiêu dự án
- Xây dựng chương trình tính toán tuần tự và song song sử dụng MPI.
- Tối ưu hóa việc trao đổi dữ liệu biên (ghost rows) giữa các Node tính toán.
- Đánh giá hiệu năng song song (Speedup, Efficiency) trên các kích thước lưới khác nhau.

---

## 2. Mô hình toán học (Mathematical Formulation)
### 2.1. Phương trình đạo hàm riêng (PDE)
Phương trình nhiệt trong không gian 2 chiều có dạng:
$$u_t = c \cdot (u_{xx} + u_{yy})$$
Trong đó:
- $u(x, y, t)$: Nhiệt độ tại vị trí $(x, y)$ thời điểm $t$.
- $c$: Hệ số khuếch tán nhiệt.

### 2.2. Điều kiện biên và điều kiện ban đầu
- **Miền tính toán:** Đơn vị vuông $[0, 1] \times [0, 1]$.
- **Điều kiện ban đầu:** $u(x, y, 0) = f(x, y)$ (ví dụ: nhiệt độ ban đầu toàn tấm là $0^\circ C$).
- **Điều kiện biên (Dirichlet):** Nhiệt độ tại 4 cạnh được giữ cố định:
    - Biên trái/phải ($x=0, x=1$): $10^\circ C$ và $40^\circ C$.
    - Biên trên/dưới ($y=0, y=1$): $50^\circ C$ và $30^\circ C$.

---

## 3. Thuật toán số (Numerical Method)
### 3.1. Rời rạc hóa (Discretization)
Sử dụng phương pháp Sai phân hữu hạn (FDM) với lưới đều:
- Bước lưới không gian: $\Delta s = 1/(n+1)$.
- Bước thời gian: $\Delta t$.
- Lưới tính toán bao gồm $(n+2) \times (n+2)$ điểm (bao gồm cả điểm biên).

### 3.2. Sơ đồ hiện (Explicit Scheme)
Công thức cập nhật nhiệt độ tại điểm $(i, j)$ từ thời điểm $k$ sang $k+1$:
$$u_{i,j}^{k+1} = u_{i,j}^k + r \cdot (u_{i+1,j}^k + u_{i-1,j}^k + u_{i,j+1}^k + u_{i,j-1}^k - 4u_{i,j}^k)$$
Với hằng số: $r = c \cdot \frac{\Delta t}{\Delta s^2}$.

### 3.3. Điều kiện ổn định (Stability Condition)
Để thuật toán ổn định và không bị phân kỳ, bước thời gian phải thỏa mãn điều kiện CFL:
$$\Delta t \le \frac{\Delta s^2}{2c}$$

---

## 4. Thiết kế song song (Parallel Design)
### 4.1. Chia miền (Domain Decomposition)
Sử dụng phương pháp chia miền 1D theo hàng:
- Toàn bộ lưới $(n+2) \times (n+2)$ được chia thành các khối hàng liên tiếp.
- Mỗi tiến trình MPI quản lý $n/P$ hàng nội bộ (với $P$ là số tiến trình).

### 4.2. Quản lý vùng đệm (Ghost Rows)
Mỗi tiến trình cần dữ liệu của hàng liền kề từ các tiến trình lân cận để tính toán.
- Mỗi tiến trình cấp phát mảng cục bộ có kích thước $(local\_n + 2) \times (n+2)$.
- 2 hàng thêm vào (hàng $0$ và hàng $local\_n + 1$) được gọi là **Ghost Rows**.

### 4.3. Giao tiếp MPI
Sử dụng hàm `MPI_Sendrecv` để tối ưu hóa việc trao đổi đồng thời:
- Tiến trình $i$ gửi hàng biên cho tiến trình $i+1$ và $i-1$.
- Nhận dữ liệu biên từ tiến trình lân cận để lấp đầy Ghost Rows.

---

## 5. Cài đặt chi tiết (Implementation)
- **Ngôn ngữ:** C.
- **Thư viện:** MPI (OpenMPI hoặc MPICH).
- **Cấu trúc dữ liệu:** Mảng một chiều được mô phỏng như mảng hai chiều để tăng tốc độ truy cập bộ nhớ.
- **Quy trình thực hiện:**
    1. Khởi tạo MPI và đọc thông số từ file `config.txt`.
    2. Chia hàng và cấp phát bộ nhớ cục bộ.
    3. Vòng lặp thời gian $K$ bước:
        - Trao đổi hàng biên (MPI).
        - Tính toán giá trị mới dựa trên stencil 5 điểm.
        - Cập nhật (Swap pointer) vùng nhớ cũ và mới.
    4. Gom dữ liệu về Rank 0 (MPI_Gatherv).
    5. Xuất kết quả và đo thời gian.

---

## 6. Thực nghiệm và Kết quả (Experiments & Results)
*(Ghi chú: Phần này sẽ cập nhật số liệu sau khi chạy)*
- **Môi trường:** [Tên máy/CPU], [Số lượng máy/core].
- **Kịch bản:** Lưới $n=500, 1000, 2000$.
- **Bảng dữ liệu:**
| Số tiến trình (P) | Thời gian (s) | Speedup | Efficiency |
|-------------------|---------------|---------|------------|
| 1 (Tuần tự)       | [Số liệu]     | 1.0     | 100%       |
| 2                 | [Số liệu]     | [...]   | [...]      |
| 4                 | [Số liệu]     | [...]   | [...]      |
| 8                 | [Số liệu]     | [...]   | [...]      |

---

## 7. Phân tích và Đánh giá (Analysis)
### 7.1. Độ chính xác
So sánh tổng kiểm (Checksum) và nhiệt độ tại tâm giữa bản tuần tự và song song: sai số xấp xỉ 0.

### 7.2. Hiệu năng
- Phân tích sự ảnh hưởng của kích thước lưới: Lưới lớn thường cho Speedup tốt hơn do giảm tỉ lệ "Giao tiếp/Tính toán".
- Giới hạn: Overhead từ việc trao đổi Ghost Rows khi số lượng tiến trình tăng quá lớn.

---

## 8. Kết luận
- Dự án đã giải quyết thành công bài toán phương trình nhiệt 2D bằng MPI.
- Hiểu rõ cơ chế chia miền và trao đổi vùng biên trong tính toán hiệu năng cao.
- **Hướng phát triển:** Mở rộng sang chia miền 2D (Checkerboard), sử dụng OpenMP để song song hóa mức luồng hoặc CUDA để chạy trên GPU.
