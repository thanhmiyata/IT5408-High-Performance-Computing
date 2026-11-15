# ĐẠI HỌC BÁCH KHOA HÀ NỘI
### TRƯỜNG CÔNG NGHỆ THÔNG TIN VÀ TRUYỀN THÔNG

---

# BÁO CÁO BÀI TẬP TUẦN 2
## HỌC PHẦN: IT5408 – TÍNH TOÁN HIỆU NĂNG CAO (HPC)
**Chủ đề:** Lập trình song song với OpenMP, Giải thuật Nhánh & Cận (TSP), Đếm số nguyên tố, Phân tích đồ thị và Nhân ma trận

* **Giảng viên hướng dẫn:** TS. Vũ Văn Thiệu, TS. Nguyễn Quốc Tuấn
* **Học viên thực hiện:** Phạm Tiến Thành
* **Mã số học viên (MSSV):** 20242422M
* **Khóa đào tạo:** Cao học Công nghệ Thông tin
* **Thời gian hoàn thành:** 15/11/2026

---

## MỤC LỤC
1. [BÀI TOÁN 1: TRAVELING SALESMAN PROBLEM (TSP)](#bài-toán-1--traveling-salesman-problem-tsp)
   * [1.1. Mô hình toán học và Thuật toán Nhánh & Cận](#11-mô-hình-toán-học-và-thuật-toán-nhánh--cận)
   * [1.2. Kỹ thuật song song hóa với OpenMP](#12-kỹ-thuật-song-song-hóa-với-openmp)
   * [1.3. Bảng kết quả thực nghiệm trên các bộ dữ liệu](#13-bảng-kết-quả-thực-nghiệm-trên-các-bộ-dữ-liệu)
2. [BÀI TOÁN 2: ĐẾM SỐ NGUYÊN TỐ ĐẾN N VÀ ĐỐI CHIẾU ĐỊNH LUẬT AMDAHL](#bài-toán-2--đếm-số-nguyên-tố-đến-n-và-đối-chiếu-định-luật-amdahl)
   * [2.1. Thuật toán và Cơ chế song song](#21-thuật-toán-và-cơ-chế-song-song)
   * [2.2. Bảng số liệu thực nghiệm Runtime và Speedup](#22-bảng-số-liệu-thực-nghiệm-runtime-và-speedup)
   * [2.3. Phân tích đối chiếu Amdahl Potential Speedup](#23-phân-tích-đối-chiếu-amdahl-potential-speedup)
3. [BÀI TOÁN 3: ĐẾM SỐ ĐƯỜNG ĐI ĐỘ DÀI K TRONG ĐỒ THỊ](#bài-toán-3--đếm-số-đường-đi-độ-dài-k-trong-đồ-thị)
   * [3.1. Cơ sở lý thuyết ma trận kề](#31-cơ-sở-lý-thuyết-ma-trận-kề)
   * [3.2. Thuật toán lũy thừa ma trận nhị phân song song](#32-thuật-toán-lũy-thừa-ma-trận-nhị-phân-song-song)
   * [3.3. Bảng kết quả thực nghiệm trên các đồ thị](#33-bảng-kết-quả-thực-nghiệm-trên-các-đồ-thị)
4. [BÀI TOÁN BỔ SUNG: NHÂN HAI MA TRẬN A(m x n) VÀ B(n x p)](#bài-toán-bổ-sung--nhân-hai-ma-trận-am-x-n-và-bn-x-p)
   * [4.1. Kiến trúc 4 phương pháp: Serial, std::thread, OpenMP, MPI](#41-kiến-trúc-4-phương-pháp-serial-stdthread-openmp-mpi)
   * [4.2. Bảng số liệu so sánh hiệu năng và năng lực tính toán (GFLOPS)](#42-bảng-số-liệu-so-sánh-hiệu-năng-và-năng-lực-tính-toán-gflops)
   * [4.3. Phân tích sự đánh đổi giữa Shared Memory và Distributed Memory](#43-phân-tích-sự-đánh-đổi-giữa-shared-memory-và-distributed-memory)
5. [KẾT LUẬN](#kết-luận)

---

# BÀI TOÁN 1 – TRAVELING SALESMAN PROBLEM (TSP)

### 1.1. Mô hình toán học và Thuật toán Nhánh & Cận
Cho đồ thị đầy đủ gồm $n$ đỉnh (thành phố) và ma trận trọng số $C = (c_{ij})$, trong đó $c_{ij} > 0$ là chi phí di chuyển từ thành phố $i$ đến thành phố $j$. Mục tiêu là tìm một chu trình $(v_1, v_2, \dots, v_n, v_1)$ đi qua mỗi đỉnh đúng một lần sao cho tổng chi phí:
$$f(v) = c_{v_n, v_1} + \sum_{i=1}^{n-1} c_{v_i, v_{i+1}} \to \min$$

Thuật toán Nhánh và Cận (Branch and Bound) tìm nghiệm tối ưu toàn cục bằng cách duyệt không gian trạng thái theo cây quyết định kết hợp kỹ thuật cắt tỉa:
* **Cận dưới (Lower Bound):** Tại mỗi nút trạng thái khi đã đi qua tập đỉnh $S \subset V$ với chi phí hiện tại là $cost_{\text{curr}}$, cận dưới của toàn bộ cây con được ước lượng bằng:
  $$\text{Bound} = cost_{\text{curr}} + \sum_{u \notin S} \min_{v \ne u} \{ c_{uv} \}$$
* **Quy tắc cắt tỉa (Pruning):** Nếu $\text{Bound} \ge best\_cost$ (kỷ lục chi phí nhỏ nhất tìm được hiện tại), toàn bộ nhánh con lập tức bị loại bỏ.

### 1.2. Kỹ thuật song song hóa với OpenMP
* **Phân rã công việc (Task Decomposition):** Cố định thành phố xuất phát là $1$ (đỉnh $0$). Nhánh tìm kiếm mức $1$ và mức $2$ (chọn đỉnh thứ hai $v_1$ và đỉnh thứ ba $v_2$) được phân phối song song cho các luồng OpenMP thông qua `#pragma omp parallel for schedule(dynamic)`.
* **Đồng bộ hóa kỷ lục (Global Upper Bound Sharing):** Biến kỷ lục $best\_cost$ được chia sẻ chung. Khi bất kỳ luồng nào tìm thấy chu trình có chi phí thấp hơn, nó cập nhật vào biến toàn cục thông qua vùng găng `#pragma omp critical`. Các luồng khác đọc được cận trên mới này và tự động cắt tỉa các nhánh kém hơn sớm hơn rất nhiều so với chạy tuần tự.

### 1.3. Bảng kết quả thực nghiệm trên các bộ dữ liệu
Thử nghiệm trên vi xử lý Apple M3 (8 nhân) với các file dữ liệu trong đề bài:

| File dữ liệu | Số đỉnh ($n$) | Chi phí nhỏ nhất ($\min$) | Chu trình tối ưu (1-based index) | Thời gian 1 luồng ($s$) | Thời gian 4 luồng ($s$) | Speedup ($4T$) |
| :--- | :---: | :---: | :--- | :---: | :---: | :---: |
| **tsp7.inp** | $4$ | **$7$** | $1 \to 2 \to 4 \to 3 \to 1$ | $0.000001\text{ s}$ | $0.000001\text{ s}$ | $1.00\text{x}$ |
| **tsp12.inp** | $12$ | **$192$** | $1 \to 4 \to 5 \to 11 \to 9 \to 12 \to 7 \to 2 \to 10 \to 3 \to 8 \to 6 \to 1$ | $0.002824\text{ s}$ | $0.001250\text{ s}$ | **$2.26\text{x}$** |
| **tsp13.inp** | $13$ | **$144$** | $1 \to 11 \to 7 \to 12 \to 8 \to 6 \to 13 \to 5 \to 9 \to 4 \to 2 \to 10 \to 3 \to 1$ | $0.003822\text{ s}$ | $0.001859\text{ s}$ | **$2.06\text{x}$** |
| **tsp14.inp** | $14$ | **$118$** | $1 \to 4 \to 6 \to 2 \to 9 \to 8 \to 7 \to 12 \to 3 \to 14 \to 5 \to 10 \to 11 \to 13 \to 1$ | $0.004457\text{ s}$ | $0.002665\text{ s}$ | **$1.67\text{x}$** |
| **tsp15.inp** | $15$ | **$166$** | $1 \to 4 \to 14 \to 10 \to 12 \to 9 \to 13 \to 7 \to 15 \to 6 \to 2 \to 8 \to 3 \to 11 \to 5 \to 1$ | $0.130907\text{ s}$ | $0.037650\text{ s}$ | **$3.48\text{x}$** |

**Nhận xét:**
* Với $n = 15$ (`tsp15.inp`), không gian trạng thái lên tới $\approx 15! / 2 \approx 6.5 \times 10^{11}$ đường đi. Nhờ kết hợp cận dưới chặt chẽ và song song hóa OpenMP đa luồng, thời gian thực thi giảm từ **$0.131\text{ s}$** xuống chỉ còn **$0.038\text{ s}$** (đạt Speedup **$3.48\text{ lần}$** trên 4 luồng).

---

# BÀI TOÁN 2 – ĐẾM SỐ NGUYÊN TỐ ĐẾN N VÀ ĐỐI CHIẾU ĐỊNH LUẬT AMDAHL

### 2.1. Thuật toán và Cơ chế song song
* Bài toán yêu cầu tính giá trị hàm đếm số nguyên tố:
  $$\pi(N) = \left|\{p \le N \mid p \text{ là số nguyên tố}\}\right|$$
* **Thuật toán tối ưu hóa:** Kiểm tra số nguyên tố bằng phương pháp Trial Division dạng $6k \pm 1$ với độ phức tạp $O(\sqrt{m})$ cho mỗi số $m$.
* **Song song hóa OpenMP:** Các số trong đoạn $[1, N]$ được kiểm tra độc lập (Embarrassingly Parallel). Sử dụng chỉ thị:
  ```cpp
  #pragma omp parallel for reduction(+:count) schedule(guided, 1000) num_threads(num_threads)
  for (long long i = 1; i <= n; ++i) { ... }
  ```
  Chiến lược `schedule(guided)` giúp tự động cân bằng tải linh hoạt vì việc kiểm tra các số lớn đòi hỏi nhiều vòng lặp hơn các số nhỏ.

### 2.2. Bảng số liệu thực nghiệm Runtime và Speedup

Khảo sát với $N$ từ $1$ triệu đến $100$ triệu trên Apple M3 (8 nhân):

| Kích thước bài toán ($N$) | Kết quả $\pi(N)$ | Số luồng ($P$) | Runtime thực tế ($s$) | Hệ số Speedup $S(P)$ | Hiệu suất song song $E(P)$ |
| :--- | :---: | :---: | :---: | :---: | :---: |
| **$1$ Triệu ($10^6$)** | **$78,498$** | $1$ (Serial) | $0.0302\text{ s}$ | $1.00\text{x}$ | $100.0\%$ |
| | | $2$ luồng | $0.0105\text{ s}$ | $2.88\text{x}$ | $144.2\%$ |
| | | $4$ luồng | $0.0048\text{ s}$ | **$6.29\text{x}$** | $157.3\%$ |
| | | $8$ luồng | $0.0040\text{ s}$ | **$7.63\text{x}$** | $95.4\%$ |
| **$10$ Triệu ($10^7$)** | **$664,579$** | $1$ (Serial) | $0.3754\text{ s}$ | $1.00\text{x}$ | $100.0\%$ |
| | | $2$ luồng | $0.1915\text{ s}$ | $1.96\text{x}$ | $98.0\%$ |
| | | $4$ luồng | $0.1047\text{ s}$ | **$3.58\text{x}$** | $89.6\%$ |
| | | $8$ luồng | $0.0860\text{ s}$ | **$4.36\text{x}$** | $54.6\%$ |
| **$50$ Triệu ($5 \times 10^7$)** | **$3,001,134$** | $1$ (Serial) | $3.5807\text{ s}$ | $1.00\text{x}$ | $100.0\%$ |
| | | $2$ luồng | $1.8326\text{ s}$ | $1.95\text{x}$ | $97.7\%$ |
| | | $4$ luồng | $1.0352\text{ s}$ | **$3.46\text{x}$** | $86.5\%$ |
| | | $8$ luồng | $0.8851\text{ s}$ | **$4.05\text{x}$** | $50.6\%$ |
| **$100$ Triệu ($10^8$)** | **$5,761,455$** | $1$ (Serial) | $9.6016\text{ s}$ | $1.00\text{x}$ | $100.0\%$ |
| | | $2$ luồng | $5.0390\text{ s}$ | **$1.91\text{x}$** | $95.3\%$ |
| | | $4$ luồng | $2.7260\text{ s}$ | **$3.52\text{x}$** | $88.1\%$ |
| | | $8$ luồng | $2.2838\text{ s}$ | **$4.20\text{x}$** | $52.6\%$ |

### 2.3. Phân tích đối chiếu Amdahl Potential Speedup
* **Ước lượng tham số:** Ở quy mô lớn $N = 100$ triệu ($T_1 = 9.6016\text{ s}, T_4 = 2.7260\text{ s}$), hệ số tăng tốc tại $4$ luồng đạt $S(4) = 3.52\text{x}$. Theo Định luật Amdahl:
  $$S(P) = \frac{1}{(1 - p) + \dfrac{p}{P}} \implies p \approx \mathbf{95.45\%}, \quad s = 1 - p \approx \mathbf{4.55\%}$$
  Giới hạn tăng tốc tối đa lý thuyết (Amdahl Potential Speedup khi $P \to \infty$) là:
  $$S_{\max} = \frac{1}{s} = \frac{1}{0.0455} \approx \mathbf{21.98\text{ lần}}$$
* **Hiện tượng tại 4 luồng và 8 luồng:**
  * Tại $P = 4$ luồng, tốc độ tăng gấp **$3.52$ lần** (hiệu suất $88.1\%$), cả 4 luồng đều hoạt động trên 4 lõi hiệu năng cao (P-cores).
  * Khi mở rộng lên $P = 8$ luồng, speedup đạt **$4.20\text{ lần}$**. Sự suy giảm hiệu suất (từ $88.1\%$ xuống $52.6\%$) phản ánh trung thực hiện tượng kiến trúc không đồng nhất (Heterogeneous Architecture): $4$ nhân E-cores chạy chậm hơn khiến các nhân P-cores phải dừng đợi đồng bộ tại cuối vòng lặp, minh chứng rõ nét cho hiện tượng thắt cổ chai trong tính toán đa lõi thực tế.

---

# BÀI TOÁN 3 – ĐẾM SỐ ĐƯỜNG ĐI ĐỘ DÀI K TRONG ĐỒ THỊ

### 3.1. Cơ sở lý thuyết ma trận kề
* Cho đồ thị $G = (V, E)$. Ma trận kề $A$ kích thước $V \times V$ có phần tử:
  $$A_{u, v} = \begin{cases} 1 & \text{nếu } (u, v) \in E \\ 0 & \text{nếu } (u, v) \notin E \end{cases}$$
* **Định lý toán học:** Số đường đi độ dài đúng $k$ từ đỉnh $u$ đến đỉnh $v$ bằng chính xác phần tử tại hàng $u$, cột $v$ của ma trận lũy thừa $A^k$, tức $(A^k)_{u, v}$.
* Do đó, tổng số đường đi độ dài $k$ trong toàn bộ đồ thị $G$ được xác định bởi:
  $$\text{TotalPaths}(k) = \sum_{u=0}^{V-1} \sum_{v=0}^{V-1} (A^k)_{u, v}$$

### 3.2. Thuật toán lũy thừa ma trận nhị phân song song
* Thay vì nhân liên tiếp $k$ lần tốn $O(k \cdot V^3)$, thuật toán áp dụng phương pháp **Lũy thừa nhị phân (Binary Exponentiation)** với độ phức tạp giảm xuống $O(\log k \cdot V^3)$.
* Mỗi phép nhân hai ma trận $V \times V$ được song song hóa theo dải hàng bằng OpenMP với thứ tự vòng lặp tối ưu Cache $I \to K \to J$:
  ```cpp
  #pragma omp parallel for schedule(static) num_threads(num_threads)
  for (int i = 0; i < V; ++i) {
      for (int k = 0; k < V; ++k) {
          unsigned long long a_ik = A[i][k];
          if (a_ik == 0) continue;
          for (int j = 0; j < V; ++j) {
              C[i][j] += a_ik * B[k][j];
          }
      }
  }
  ```

### 3.3. Bảng kết quả thực nghiệm trên các đồ thị

| File đồ thị | Số đỉnh ($V$) | Số cạnh ($E$) | Độ dài đường đi ($k$) | Tổng số đường đi độ dài $k$ | Thời gian tính toán ($s$) |
| :--- | :---: | :---: | :---: | :---: | :---: |
| **graph01.txt** | $6$ | $20$ | $k = 3$ | **$224$** | $0.000068\text{ s}$ |
| **graph02.txt** | $6$ | $30$ | $k = 3$ | **$750$** | $0.000052\text{ s}$ |
| **graph03.txt** | $13$ | $156$ | $k = 4$ | **$269,568$** | $0.000053\text{ s}$ |
| **graph04.txt** | $14$ | $182$ | $k = 4$ | **$399,854$** | $0.000049\text{ s}$ |
| **graph05.txt** | $15$ | $210$ | $k = 5$ | **$8,067,360$** | $0.000065\text{ s}$ |

---

# BÀI TOÁN BỔ SUNG – NHÂN HAI MA TRẬN A(m x n) VÀ B(n x p)

Theo yêu cầu từ file `baitap_additional.docx`: *"Nhân 2 ma trận A (m x n) và B (n x p) bằng: tuần tự, song song với C++ thread, song song với openmp, song song với mpi. Đo hiệu năng và so sánh các chương trình."*

### 4.1. Kiến trúc 4 phương pháp: Serial, std::thread, OpenMP, MPI
1. **Tuần tự (Serial):** Chuẩn hóa theo thứ tự duyệt $I \to K \to J$ (Cache-friendly) để tận dụng dòng Cache L1/L2 và bộ vector hóa tự động của trình biên dịch.
2. **C++ `std::thread`:** Chia đều dải hàng của ma trận $A$ cho $T$ luồng, mỗi luồng tính toán độc lập vùng hàng tương ứng trong $C$, không xảy ra xung đột ghi (No Data Race).
3. **OpenMP:** Sử dụng `#pragma omp parallel for schedule(static)` trên vòng lặp ngoài cùng $i$, tận dụng thread pool có sẵn giúp triệt tiêu chi phí tạo luồng.
4. **MPI (Distributed Memory):** Rank 0 broadcast toàn bộ $B$ (`MPI_Bcast`), phân phối các dải hàng của $A$ qua `MPI_Scatterv`, mỗi rank tính toán cục bộ phần ma trận $local\_C$, sau đó thu gom kết quả về Rank 0 qua `MPI_Gatherv`.

### 4.2. Bảng số liệu so sánh hiệu năng và năng lực tính toán (GFLOPS)

Thực nghiệm đo đạc với ma trận cỡ vừa $512 \times 512$ và cỡ lớn $1024 \times 1024$:

| Kích thước ($m \times n \times p$) | Phương pháp | Số luồng / Rank | Thời gian ($s$) | Speedup | Năng lực tính toán (GFLOPS) | Checksum kiểm tra |
| :--- | :--- | :---: | :---: | :---: | :---: | :---: |
| **$512 \times 512 \times 512$** | **Tuần tự (Serial)** | $1$ | $0.0386\text{ s}$ | $1.00\text{x}$ | $6.95\text{ GFLOPS}$ | $3.3115 \times 10^7$ |
| | C++ `std::thread` | $4$ | $0.0064\text{ s}$ | **$6.05\text{x}$** | $41.87\text{ GFLOPS}$ | $3.3115 \times 10^7$ |
| | OpenMP | $4$ | $0.0066\text{ s}$ | **$5.88\text{x}$** | $40.88\text{ GFLOPS}$ | $3.3115 \times 10^7$ |
| | OpenMPI | $4$ | $0.0202\text{ s}$ | **$1.92\text{x}$** | $13.29\text{ GFLOPS}$ | $3.3115 \times 10^7$ |
| **$1024 \times 1024 \times 1024$** | **Tuần tự (Serial)** | $1$ | $0.1504\text{ s}$ | $1.00\text{x}$ | $14.28\text{ GFLOPS}$ | $2.2498 \times 10^8$ |
| | C++ `std::thread` | $4$ | $0.0445\text{ s}$ | **$3.38\text{x}$** | $48.26\text{ GFLOPS}$ | $2.2498 \times 10^8$ |
| | OpenMP | $4$ | $0.0431\text{ s}$ | **$3.49\text{x}$** | **$49.83\text{ GFLOPS}$** | $2.2498 \times 10^8$ |
| | OpenMPI | $4$ | $0.1088\text{ s}$ | **$1.38\text{x}$** | $19.74\text{ GFLOPS}$ | $2.2498 \times 10^8$ |

*(Ghi chú: Giá trị Checksum của cả 4 phương pháp đều trùng khớp tuyệt đối, khẳng định tính đúng đắn toán học).*

### 4.3. Phân tích sự đánh đổi giữa Shared Memory và Distributed Memory
* **OpenMP và `std::thread`:** Đạt hiệu năng tối ưu trên hệ thống đa lõi dùng chung bộ nhớ (Shared Memory). Do toàn bộ dữ liệu nằm trên RAM chung, các luồng chỉ việc đọc dữ liệu mà không tốn chi phí sao chép, năng lực tính toán đạt xấp xỉ **$50\text{ GFLOPS}$**.
* **MPI:** Trên một máy tính cá nhân đơn lẻ, MPI bị suy giảm tốc độ do chi phí truyền thông bộ nhớ (IPC Memory Copy Overhead) của các hàm tập thể `MPI_Bcast` và `MPI_Scatterv`/`MPI_Gatherv`. Tuy nhiên, MPI là giải pháp duy nhất có khả năng mở rộng trên các hệ thống cụm siêu máy tính (Clusters) gồm hàng nghìn nút mạng độc lập khi kích thước ma trận vượt quá dung lượng RAM của một máy đơn.

---

# KẾT LUẬN

1. Đã giải quyết trọn vẹn cả 3 bài toán chính trong `baitap.pdf` (TSP, Đếm số nguyên tố, Đếm số đường đi đồ thị) và bài toán bổ sung trong `baitap_additional.docx` (Nhân 2 ma trận).
2. Kiểm chứng thành công thuật toán Nhánh & Cận song song cho bài toán TSP, giải chính xác các trường hợp $n \le 15$ trong thời gian mili-giây.
3. Đo lường thực nghiệm Runtime đếm số nguyên tố quy mô lớn ($N = 10^8$), qua đó đối chiếu và giải thích sâu sắc sự phù hợp cũng như giới hạn thực tế của **Định luật Amdahl** trên vi kiến trúc đa lõi hiện đại.
4. Đánh giá toàn diện sự khác biệt cốt lõi giữa hai mô hình lập trình HPC kinh điển: **Shared Memory** (OpenMP, std::thread) và **Distributed Memory** (MPI).
