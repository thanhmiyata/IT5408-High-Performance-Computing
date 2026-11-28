# KỊCH BẢN THUYẾT TRÌNH CHI TIẾT (CHUẨN 15 PHÚT)
**Đề tài:** Giải thuật và Cài đặt Chương trình Song song Bài toán có Sự phụ thuộc Dữ liệu — Lattice Boltzmann Method (D2Q9)  
**Nhóm thực hiện:** Nhóm 1 — IT5408 Tính toán Hiệu năng Cao

> **Lưu ý nhanh trước khi lên bục (từ rà soát Slide 1–20):**
> - **Slide 15 & 16:** Slide hiện có bảng số tham khảo cũ và 2 chỗ ghi `[Velocity Field Plot]` / `[Performance Charts]` chưa chèn ảnh. Đến đoạn này bạn dẫn mượt sang **Demo trực tiếp trên Terminal + Web UI (`lbm_demo_ui.html`)** để thầy xem số chạy thật và biểu đồ thật luôn!
> - **Slide 10 & 14:** Slide minh họa theo sơ đồ tổng quát (push scheme / điều kiện `rank > 0`), còn trong code thực tế chạy ở file `src/lbm_mpi.c` nhóm đã cải tiến sang pull scheme + vòng tuần hoàn `(rank ± 1 + size) % size` để biên X tuần hoàn khớp tuyệt đối bản tuần tự. Nếu thầy có hỏi thì nói rõ điểm cải tiến này.

---

## CHẶNG 1: MỞ ĐẦU & GIỚI THIỆU (Phút 0:00 – 1:30)

### Slide 1 – 2 (Trang bìa & Tiêu đề) — *~20 giây*
> "Em xin kính chào Thầy và các anh/chị học viên. Em là **[Tên bạn]**, đại diện cho **Nhóm 1** lớp IT5408 – Tính toán hiệu năng cao.  
> Hôm nay nhóm em xin trình bày đồ án cuối kỳ với đề tài: **'Giải thuật và cài đặt chương trình song song bài toán có sự phụ thuộc dữ liệu – Ứng dụng phương pháp Lattice Boltzmann (D2Q9) bằng MPI'**."

### Slide 3 (Thành viên Nhóm 1) — *~25 giây*
> "Nhóm em gồm 6 thành viên với phân công rõ ràng:  
> - Bạn Phạm Tiến Thành làm trưởng nhóm, phụ trách cài đặt song song MPI và benchmark;  
> - Bạn Trần Thị Hạnh phụ trách mô hình toán D2Q9 và điều kiện biên;  
> - Bạn Bùi Minh Đức cài đặt bản tuần tự C99 và kiểm thử checksum;  
> - Bạn Nguyễn Thị Huê xử lý dữ liệu và trực quan hóa;  
> - Bạn Nguyễn Đình Hiếu thiết kế slide và phân tích scaling;  
> - Bạn Đỗ Minh Hiếu phụ trách báo cáo kỹ thuật và phân tích overhead."

### Slide 4 (Nội dung trình bày) — *~20 giây*
> "Bài trình bày của nhóm gồm 8 phần chính: từ mô hình toán học, phân tích sự phụ thuộc dữ liệu, thiết kế giải thuật song song SPMD và truyền thông MPI, cho đến phần demo chạy trực tiếp và phân tích kết quả thực nghiệm đa kích thước lưới."

---

## CHẶNG 2: MÔ HÌNH TOÁN HỌC & SỰ PHỤ THUỘC DỮ LIỆU (Phút 1:30 – 4:30)

### Slide 5 (Phương trình Lattice Boltzmann – LBE) — *~50 giây*
> "Trước hết là phần 1: Mô hình toán học. Khác với các phương pháp truyền thống như sai phân hữu hạn (FDM) hay phần tử hữu hạn (FEM) giải trực tiếp hệ phương trình vi phân Navier-Stokes phi tuyến, **Lattice Boltzmann Method (LBM)** xuất phát từ lý thuyết động học phân tử: thay vì tính trực tiếp vận tốc hay áp suất, ta theo dõi **hàm phân bố mật độ hạt $f_i(x, t)$**.  
> Phương trình LBE với xấp xỉ va chạm BGK tách gọn thành **2 bước lặp** ở mỗi bước thời gian:  
> 1. **Collision (Va chạm):** tính giá trị trung gian $f_i^*$ hướng về trạng thái cân bằng $f_i^{eq}$ với hệ số thư giãn $\omega$.  
> 2. **Streaming (Lan truyền):** dịch chuyển các hạt $f_i^*$ sang các nút lưới lân cận theo vector vận tốc $c_i$."

### Slide 6 (Mô hình D2Q9 Stencil) — *~45 giây*
> "Nhóm em sử dụng mô hình kinh điển **D2Q9**, tức không gian **2 chiều (D2)** với **9 hướng vận tốc rời rạc (Q9)** tại mỗi điểm lưới như trên sơ đồ:  
> - Hướng `0` đứng yên ở tâm với trọng số $w_0 = 4/9$;  
> - 4 hướng chính (lên, xuống, trái, phải) có trọng số $w_{1..4} = 1/9$;  
> - 4 hướng chéo có trọng số $w_{5..8} = 1/36$. Tổng trọng số bằng $1$."

### Slide 7 (Phân tích sự phụ thuộc dữ liệu) — *~50 giây*
> "Sang phần 2 – đây là trọng tâm bài toán HPC: **Phân tích sự phụ thuộc dữ liệu**.  
> - Ở bước **Collision**, phép tính tại ô $(x, y)$ chỉ dùng dữ liệu của chính ô đó → **hoàn toàn độc lập**, không có phụ thuộc dữ liệu giữa các nút.  
> - Sự phụ thuộc dữ liệu chỉ xảy ra ở bước **Streaming**: để cập nhật giá trị tại điểm $x$ ở bước $t+1$, ta cần lấy giá trị $f_i^*$ từ **8 ô láng giềng xung quanh**.  
> - Trong bản tuần tự, biên trên/dưới áp dụng tường cứng **Bounce-back** và biên trái/phải theo điều kiện tuần hoàn **Periodic**. Còn khi song song hóa trên bộ nhớ phân tán MPI, mỗi CPU giữ một phần lưới nên bắt buộc phải tạo các **Ghost Columns (cột biên ảo)** và trao đổi dữ liệu biên giữa các tiến trình ngay trước bước Streaming."

### Slide 8 (Đại lượng vĩ mô & Ký hiệu) — *~35 giây*
> "Từ 9 hàm phân bố $f_i$ tại mỗi điểm lưới, ta suy ra trực tiếp hai đại lượng vĩ mô: **mật độ $\rho$** bằng tổng các $f_i$, và **vận tốc $\mathbf{u}$** bằng tổng động lượng $f_i \mathbf{c}_i$ chia cho $\rho$. Bên phải là công thức khai triển hàm cân bằng Maxwell-Boltzmann rời rạc $f_i^{eq}$ để đưa vào bước Collision."

---

## CHẶNG 3: CÀI ĐẶT & GIẢI THUẬT SONG SONG MPI (Phút 4:30 – 8:00)

### Slide 9 & 10 (Code Collision & Streaming với Ghost Columns) — *~50 giây*
> "Về cài đặt thực tế bằng ngôn ngữ C thuần không dùng thư viện ngoài:  
> - Ở **Slide 9**, hàm `CollisionCucBo` duyệt trên từng cột cục bộ `nx_local`, tính $f^{eq}$ và cập nhật $f_{new}$ hoàn toàn tại chỗ mà không cần gọi bất kỳ hàm MPI nào.  
> - Ở **Slide 10**, hàm `StreamingCucBo` xử lý 3 trường hợp: nếu hạt chạm thành tường trên/dưới (`y < 0` hoặc `y >= ny`) thì bật ngược lại theo `opposite[i]`; nếu vượt ra ngoài ranh giới trái/phải của rank hiện tại thì đọc từ vùng đệm `ghost_left` hoặc `ghost_right` đã nhận về (lưu ý lấy đúng hàng `ys = y - cy[i]` cho cả các hướng chéo); còn lại thì lan truyền bình thường trong miền."

### Slide 11 (Chiến lược Song song hóa SPMD) — *~50 giây*
> "Nhóm áp dụng mô hình **SPMD (Single Program, Multiple Data)** với chiến lược **chia miền 1 chiều (1D Domain Decomposition) theo trục X (theo cột)**.  
> Vì sao lại cắt theo cột mà không cắt theo hàng? Vì trong C, mảng cấp phát theo thứ tự `(x * ny + y) * Q`, nên toàn bộ một **cột `x` nằm liên tiếp nhau trên bộ nhớ RAM**. Nhờ vậy khi gửi cột biên qua MPI, ta gửi thẳng một vùng nhớ liên tục kích thước `ny × 9 × 8 bytes` mà không tốn chi phí copy hay đóng gói dữ liệu.  
> Ví dụ với lưới `256 × 64` chạy trên 4 tiến trình, mỗi rank phụ trách 64 cột và chỉ cần trao đổi cột đệm ~4.5 KB với hàng xóm trái/phải."

### Slide 12 (5 Bước chính trong giải thuật) — *~35 giây*
> "Tổng thể chương trình MPI gồm 5 bước: (1) Khởi tạo bộ nhớ và trạng thái tĩnh $\rho=1, u=0$; (2) Chia miền `nx_local` và xử lý phần dư `nx % size` để cân bằng tải; (3) Phát tán cấu hình từ Rank 0 qua `MPI_Bcast`; (4) Chạy vòng lặp `10,000` bước gồm *Tính Macro → Áp biên → Collision → Trao đổi Ghost → Streaming*; và (5) Thu thập Checksum bằng `MPI_Reduce` và gom dữ liệu bằng `MPI_Gatherv` về Rank 0 xuất file."

### Slide 13 & 14 (Chi tiết truyền thông Non-blocking: `TraoDoiGhost`) — *~55 giây*
> "Điểm then chốt về truyền thông nằm ở hàm `TraoDoiGhost` (Slide 14).  
> Thay vì dùng `MPI_Send` / `MPI_Recv` dạng chặn (blocking) vốn rất dễ gây **Deadlock** khi tất cả các rank cùng gửi sang hàng xóm và chờ nhau, nhóm sử dụng cặp hàm không chặn **`MPI_Irecv` và `MPI_Isend`**, kết hợp với **`MPI_Waitall`** để đồng bộ trước khi vào bước Streaming.  
> Nhóm cũng phân tách 2 tag `100` và `101` cho luồng truyền sang phải và sang trái để tránh nhận nhầm gói tin khi chạy 2 tiến trình, và dùng chỉ số vòng `(rank ± 1 + size) % size` để khớp với biên tuần hoàn theo trục X."

---

## CHẶNG 4: THỰC NGHIỆM, LIVE DEMO & KẾT LUẬN (Phút 8:00 – 13:30)

### Slide 15 & 16 → Chuyển sang DEMO TRỰC TIẾP — *~4 phút 30 giây*
*(Nói ngắn gọn ở Slide 15–16 rồi bấm Alt+Tab sang Terminal và Web UI)*

> **Dẫn trên Slide 15–16 (~20 giây):**  
> "Trên Slide 15–16 là bảng tổng hợp và các chỉ số đánh giá chính gồm **MLUPS** (số triệu điểm lưới cập nhật mỗi giây), **Speedup**, **Efficiency** và đối soát **Checksum**. Thay vì chỉ nhìn số tĩnh trên slide, em xin phép chạy **demo trực tiếp trên máy thật** và mở giao diện trực quan hóa luồng nước để Thầy quan sát kết quả mới nhất của nhóm."

> **1. Chạy Live trên Terminal (~1 phút):**  
> *(Gõ lệnh: `./lbm_serial -c config.txt && mpirun -np 4 ./lbm_mpi -c config.txt`)*  
> "Ở đây em chạy cùng cấu hình lưới `256 × 64`, `10,000` bước trên chip 8 nhân:  
> - Bản tuần tự `lbm_serial` chạy mất khoảng **1.58 giây**, đạt **~103 MLUPS**.  
> - Bản song song `lbm_mpi` với 4 tiến trình chỉ mất khoảng **0.51 giây**, đạt **~323 MLUPS** — tức tăng tốc **3.12 lần**, hiệu suất đạt **78.1%**.  
> - Và Thầy có thể thấy ngay phần **Checksum**: tổng mật độ `Sum(rho)` của Serial là `16383.9999999909`, còn của MPI là `16383.9999999910`, `Sum(ux)` là `1527.0476269322` — **khớp nhau tuyệt đối đến 10 chữ số thập phân**, chứng minh thuật toán song song bảo toàn hoàn toàn khối lượng và kết quả vật lý."

> **2. Mở giao diện UI `lbm_demo_ui.html` (~2 phút):**  
> *(Bấm `Về Step 0` → `Play`)*  
> "Để quan sát dòng nước thực sự chảy như thế nào trong quá trình tính toán, nhóm em đã thêm tham số xuất chuỗi snapshot theo chu kỳ và xây dựng giao diện này:  
> - Ở **Step 0**, nước trong ống đứng yên hoàn toàn, chỉ có hai biên đầu vào/ra đẩy vận tốc $u_0 = 0.1$.  
> - Khi cho chạy qua các step `1,000 → 5,000`, xung lượng lan truyền vào giữa kênh. Nhờ điều kiện *no-slip* ở tường trên và dưới, dòng sát tường bị giữ bằng 0, còn ở tâm ống chảy nhanh dần đạt `~0.141`, hình thành đúng **mặt cắt Parabol Poiseuille** khớp với đường nét đứt lý thuyết màu đỏ ở biểu đồ bên phải.  
> - Thầy có thể thấy **4 vạch chia `Rank 0 | Rank 1 | Rank 2 | Rank 3`**: các hạt nước và trường màu chảy qua ranh giới giữa các Rank **hoàn toàn liền mạch**, không có đường gãy hay sai số biên.  
> *(Bấm chuyển nhanh sang tab `512 × 128` và `1024 × 256`)*  
> - Ngoài lưới chuẩn `256 × 64`, nhóm cũng đã mở rộng chạy kiểm chứng trên lưới lớn **`512 × 128`** (65 nghìn nút) và **`1024 × 256`** (262 nghìn nút) như hiển thị trên màn hình."

> **3. Kéo xuống phần Biểu đồ so sánh 3 lưới ngay trên UI (~1 phút 30 giây):**  
> "Nhìn vào biểu đồ Strong Scaling của 3 kích thước lưới ở phía dưới, nhóm rút ra 3 nhận xét quan trọng về kiến trúc HPC:  
> 1. Với lưới `256 × 64` (tốn `~2.36 MB`), dữ liệu nằm gọn trong **CPU Cache** nên chạy 4 tiến trình tăng tốc rất tốt (**3.12×**).  
> 2. Khi tăng lên `512 × 128` (`~9.44 MB`) và `1024 × 256` (`~37.75 MB`), Speedup tại `P = 4` giảm dần xuống **2.22×** và **1.59×**. Nguyên nhân là vì LBM thuộc nhóm bài toán **Memory-Bound** (đọc/ghi 18 biến `double` ở mỗi nút mỗi bước); khi dữ liệu tràn khỏi Cache xuống RAM dùng chung trên một máy đơn thì băng thông RAM bị bão hòa — đây chính là minh chứng thực tế vì sao bài toán LBM quy mô lớn cần chạy MPI trên cụm **Cluster nhiều node**.  
> 3. Tại `P = 8`, hiệu năng giảm do khối lượng tính toán trên mỗi cột nhỏ đi trong khi vẫn phải đồng bộ `MPI_Waitall` ở mỗi bước, cộng thêm hiệu ứng **Straggler** trên kiến trúc nhân lai (P-core phải chờ E-core)."

### Slide 17 – 20 (Kết luận, Hướng phát triển & Cảm ơn) — *~1 phút*
> *(Quay lại Slide 17 → 20)*  
> "Tóm lại, nhóm đã hoàn thành cài đặt độc lập LBM D2Q9 tuần tự và song song MPI, đảm bảo tính đúng đắn tuyệt đối và trực quan hóa đầy đủ.  
> Hạn chế của LBM là tiêu tốn bộ nhớ (9 hướng mỗi nút) và dễ gặp bottleneck băng thông RAM trên lưới lớn. Trong tương lai, hướng phát triển tự nhiên nhất là chuyển sang **chia miền 2D**, kết hợp **Hybrid MPI + OpenMP**, hoặc đưa hai bước Collision/Streaming lên **GPU (CUDA)** để đạt tăng tốc hàng chục đến hàng trăm lần.  
> Phần trình bày và demo của Nhóm 1 đến đây là kết thúc. Em xin trân trọng cảm ơn Thầy và các bạn đã lắng nghe, nhóm em rất mong nhận được các câu hỏi và góp ý từ Thầy ạ!"

---

## 3. GIẢI THÍCH CÁC CHỈ SỐ TRÊN GIAO DIỆN DEMO (`lbm_demo_ui.html`)

### A. Thông số cấu hình & Trạng thái mô phỏng (Nửa trên UI)
| Chỉ số / Thành phần | Ý nghĩa & Giá trị thực tế | Cách giải thích ngắn khi Thầy hỏi |
|:---|:---|:---|
| **`Lưới: 256×64 / 512×128 / 1024×256`** | Kích thước không gian 2D ($N_x \times N_y$) | Tương ứng `16,384`, `65,536` và `262,144` điểm lưới cần cập nhật mỗi bước. |
| **`MPI Ranks: 4`** | Số tiến trình song song MPI ($P = 4$) | Mỗi tiến trình xử lý 1/4 số cột của kênh và trao đổi qua Ghost Columns. |
| **`ω = 1.0`, `u₀ = 0.1`** | Tham số thư giãn BGK & Vận tốc đầu vào | $\omega = 1/\tau \implies$ độ nhớt $\nu = \frac{1/\omega - 0.5}{3} = 1/6 \approx 0.1667$; chọn $u_0 \le 0.1$ để số Mach nhỏ ($Ma \ll 0.3$), đảm bảo ổn định số. |
| **`Vận tốc lớn nhất \|u\|`** | $\max \sqrt{u_x^2 + u_y^2}$ trên toàn kênh | Đạt `~0.1408` (lớn hơn $u_0 = 0.1$ vì thành trên/dưới hãm dòng sát tường về 0 nên tâm ống chảy nhanh hơn để bảo toàn lưu lượng). |
| **`Vận tốc tâm ống`** | Giá trị $u_x$ tại điểm giữa $(N_x/2, N_y/2)$ | Tăng từ `0.0000` (Step 0) → `~0.1382` (Step 5,000+) rồi giữ ổn định khi dòng hội tụ. |
| **`Mật độ TB (ρ)`** | Trung bình mật độ $\rho$ toàn kênh | Luôn giữ đúng **`1.00000`** qua mọi step → chứng minh **bảo toàn khối lượng tuyệt đối**. |
| **`Trạng thái dòng`** | Pha phát triển của chất lưu | Tự động phân loại: *Khởi tạo tĩnh* (Step 0) → *Đang lan truyền* → *Parabol ổn định*. |
| **`Mặt cắt Parabol u_x(y)`** | Đồ thị vận tốc theo chiều cao $y$ tại $x = N_x/2$ | **Đường xanh** (LBM) chồng khít lên **đường đỏ nét đứt** (nghiệm giải tích Poiseuille $u_{\text{tâm}}[1 - (\frac{y - H/2}{H/2})^2]$). |
| **`Trường hiển thị` & `Ranh giới Ranks`** | Xem `\|u\|`, `u_x`, hoặc `ρ` (áp suất $p = \rho/3$) | Các vạch nét đứt `Rank 0..3` minh họa chia miền 1D; hạt nước & màu chảy qua ranh giới **liền mạch hoàn toàn**. |

### B. Chỉ số Đo hiệu năng & So sánh Đa lưới (Nửa dưới UI)
| Chỉ số | Công thức / Đơn vị | Ý nghĩa phân tích |
|:---|:---|:---|
| **`Dung lượng bộ nhớ (f+f_new)`** | $2 \times N_x \times N_y \times 9 \times 8$ bytes | `256×64` (`2.36 MB`) nằm gọn trong **CPU Cache**; `512×128` (`9.44 MB`) & `1024×256` (`37.75 MB`) tràn xuống **RAM/DRAM** gây nghẽn băng thông. |
| **`Thời gian (s)`** | Giây (`MPI_Wtime` / `clock`) | Thời gian thuần túy chạy vòng lặp chính (không tính đọc/ghi file đầu cuối). |
| **`MLUPS`** | $\frac{N_x \times N_y \times N_{\text{steps}}}{T \times 10^6}$ | Số **triệu nút lưới cập nhật mỗi giây** — thước đo chuẩn của LBM (đỉnh **`323.26 MLUPS`** tại `256×64`, $P=4$). |
| **`Speedup`** | $S_P = T_1 / T_P$ | Hệ số tăng tốc so với tuần tự 1 nhân (đạt **`3.12×`** tại `256×64`, $P=4$; đường nét đứt xám là *Lý tưởng* $S_P = P$). |
| **`Efficiency (%)`** | $E_P = \frac{S_P}{P} \times 100\%$ | Mức tận dụng mỗi CPU (đạt **`78.1%`** tại $P=4$; phần hụt ~22% là chi phí giao tiếp và đồng bộ `MPI_Waitall`). |
| **`Checksum Sum(ρ)`** | Tổng $\sum \rho$ cuối mô phỏng (`MPI_Reduce`) | Sai lệch giữa MPI và Serial `< 10⁻¹²%`, khẳng định tính đúng đắn số học tuyệt đối. |
