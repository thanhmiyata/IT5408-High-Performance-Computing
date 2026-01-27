# THIẾT KẾ SONG SONG MPI CHO LBM

## 1. Ý tưởng chia miền

### 1.1. Tại sao LBM phù hợp với song song hóa?
- **Tính cục bộ cao:** Mỗi điểm lưới chỉ cần dữ liệu từ các điểm lân cận trực tiếp.
- **Cấu trúc đều đặn:** Lưới Cartesian đều đặn, dễ phân chia.
- **Tính toán đồng nhất:** Mọi điểm lưới thực hiện cùng một phép tính.

### 1.2. Chiến lược chia miền
Sử dụng **Domain Decomposition 1D theo chiều X** (chia theo cột):
- Toàn bộ lưới `Nx × Ny` được chia thành `P` khối theo chiều X.
- Mỗi tiến trình MPI quản lý `Nx/P` cột.
- Chiều Y không chia (mỗi tiến trình có đầy đủ chiều cao).

**Ví dụ:** Lưới 256×64 với 4 tiến trình:
```
Rank 0: cột 0-63   (64 cột)
Rank 1: cột 64-127 (64 cột)
Rank 2: cột 128-191 (64 cột)
Rank 3: cột 192-255 (64 cột)
```

## 2. Vùng đệm ảo (Ghost Columns)

### 2.1. Vấn đề
Trong bước streaming, hạt di chuyển sang các điểm lân cận:
- Hạt ở cột biên phải của Rank i cần đi sang cột trái của Rank i+1.
- Hạt ở cột biên trái của Rank i cần đi sang cột phải của Rank i-1.

### 2.2. Giải pháp
Mỗi tiến trình cấp phát thêm 2 cột ảo (ghost columns):
- **Ghost Left:** Chứa dữ liệu từ Rank i-1.
- **Ghost Right:** Chứa dữ liệu từ Rank i+1.

### 2.3. Cấu trúc bộ nhớ
```
Rank i:
  Cột ảo trái | Cột 0 | Cột 1 | ... | Cột (nx_local-1) | Cột ảo phải
  [ghost_left]  [-------- nx_local cột thực ---------]  [ghost_right]
```

Kích thước mảng cục bộ:
- Không có ghost: `nx_local × ny × 9` (9 hướng)
- Có ghost: Lưu riêng 2 mảng `ghost_left[ny × 9]` và `ghost_right[ny × 9]`

## 3. Giao tiếp MPI

### 3.1. Quy trình trao đổi biên
Sau bước **Collision**, trước bước **Streaming**:

```
1. Rank i gửi cột biên phải (cột nx_local-1) cho Rank i+1
2. Rank i gửi cột biên trái (cột 0) cho Rank i-1
3. Rank i nhận ghost_left từ Rank i-1
4. Rank i nhận ghost_right từ Rank i+1
```

### 3.2. Sử dụng MPI_Irecv/Isend (Non-blocking)
Để tránh deadlock và tối ưu hiệu năng:

```c
MPI_Request reqs[4];
int rcount = 0;

// Gửi/nhận với rank trái
if (rank > 0) {
  MPI_Irecv(ghost_left, ny*Q, MPI_DOUBLE, rank-1, 100, comm, &reqs[rcount++]);
  MPI_Isend(f_new + 0*ny*Q, ny*Q, MPI_DOUBLE, rank-1, 101, comm, &reqs[rcount++]);
}

// Gửi/nhận với rank phải
if (rank < size-1) {
  MPI_Irecv(ghost_right, ny*Q, MPI_DOUBLE, rank+1, 101, comm, &reqs[rcount++]);
  MPI_Isend(f_new + (nx_local-1)*ny*Q, ny*Q, MPI_DOUBLE, rank+1, 100, comm, &reqs[rcount++]);
}

MPI_Waitall(rcount, reqs, MPI_STATUSES_IGNORE);
```

### 3.3. Xử lý biên chu kỳ (Periodic)
Nếu sử dụng điều kiện biên chu kỳ theo X:
- Rank 0 gửi cho Rank (size-1) và ngược lại.
- Hoặc xử lý trực tiếp trong hàm streaming.

## 4. Thuật toán song song

### 4.1. Khởi tạo
```
1. MPI_Init()
2. Xác định rank và size
3. Tính nx_local = nx / size (có xử lý phần dư)
4. Cấp phát bộ nhớ cục bộ:
   - f[nx_local × ny × 9]
   - f_new[nx_local × ny × 9]
   - rho[nx_local × ny]
   - ux, uy[nx_local × ny]
   - ghost_left[ny × 9], ghost_right[ny × 9]
5. Khởi tạo fᵢ = fᵢᵉᑫ với ρ=1, u=0
```

### 4.2. Vòng lặp thời gian
```
for (step = 0; step < nsteps; step++) {
  // 1. Tính các đại lượng vĩ mô
  TinhMacroCucBo(f, rho, ux, uy, nx_local, ny);
  
  // 2. Collision
  CollisionCucBo(f, f_new, rho, ux, uy, nx_local, ny, omega);
  
  // 3. Trao đổi biên
  TraoDoiGhost(f_new, nx_local, ny, ghost_left, ghost_right, rank, size);
  
  // 4. Streaming
  StreamingCucBo(f_new, f, nx_local, ny, ghost_left, ghost_right, rank, size);
  
  // 5. Áp dụng điều kiện biên (bounce-back, periodic, etc.)
}
```

### 4.3. Thu thập kết quả
Sau khi mô phỏng xong, có thể:
- **Tính toán phân tán:** Mỗi rank tính chỉ số cục bộ, sau đó dùng `MPI_Reduce` để tổng hợp.
- **Thu thập toàn bộ:** Dùng `MPI_Gather` để rank 0 thu thập toàn bộ dữ liệu và ghi file.

## 5. Phân tích hiệu năng

### 5.1. Lượng dữ liệu trao đổi
Mỗi bước thời gian, mỗi tiến trình gửi/nhận:
- **Dữ liệu gửi:** 2 cột × ny × 9 × sizeof(double) bytes
- **Ví dụ:** ny=64, 9 hướng, double=8 bytes → 2 × 64 × 9 × 8 = 9216 bytes ≈ 9 KB

### 5.2. Tỉ lệ Tính toán / Giao tiếp
```
Computation: nx_local × ny × 9 × (số phép toán)
Communication: 2 × ny × 9 × (latency + bandwidth)
```

Khi `nx_local` lớn, tỉ lệ này tốt → Hiệu năng cao.
Khi `nx_local` nhỏ (quá nhiều tiến trình), overhead giao tiếp tăng.

### 5.3. Speedup lý tưởng
```
Speedup = T_serial / T_parallel
```

Trong thực tế, speedup < P do:
- Overhead giao tiếp MPI.
- Overhead khởi tạo và đồng bộ.
- Load imbalance (nếu nx không chia hết cho P).

## 6. Tối ưu hóa

### 6.1. Overlapping Communication và Computation
Sử dụng MPI non-blocking để tính toán các điểm bên trong trong khi chờ giao tiếp biên.

### 6.2. Chia miền 2D
Thay vì chia 1D, có thể chia lưới thành các khối 2D (ví dụ: 4×4 = 16 khối).
- **Ưu điểm:** Giảm lượng dữ liệu trao đổi (tỉ lệ diện tích/chu vi tốt hơn).
- **Nhược điểm:** Code phức tạp hơn.

### 6.3. Hybrid MPI + OpenMP
Kết hợp MPI (giữa các node) và OpenMP (trong mỗi node) để tận dụng đa nhân.

## 7. Xử lý điều kiện biên trong MPI

### 7.1. Bounce-back (Tường cứng)
Xử lý cục bộ, không cần giao tiếp:
```c
if (yn < 0 || yn >= ny) {
  int opp = opposite[i];
  f[x][y][opp] += f_new[x][y][i];
}
```

### 7.2. Periodic (Chu kỳ)
- **Theo Y:** Xử lý cục bộ.
- **Theo X:** Cần giao tiếp giữa Rank 0 và Rank (size-1).

### 7.3. Zou-He (Vận tốc/Áp suất)
Chỉ rank chứa biên đầu vào/ra mới cần xử lý.

## 8. Ví dụ phân chia

### Lưới 256×64 với 4 tiến trình:

| Rank | Cột quản lý | nx_local | Số điểm | Dữ liệu (MB) |
|------|-------------|----------|---------|--------------|
| 0    | 0-63        | 64       | 4096    | 0.28         |
| 1    | 64-127      | 64       | 4096    | 0.28         |
| 2    | 128-191     | 64       | 4096    | 0.28         |
| 3    | 192-255     | 64       | 4096    | 0.28         |

**Tổng:** 16384 điểm, 1.13 MB (với 9 hướng, double precision)

## 9. Checklist triển khai

- ✅ Chia miền đều (xử lý phần dư nếu nx % size != 0)
- ✅ Cấp phát ghost columns
- ✅ Trao đổi biên bằng MPI_Irecv/Isend
- ✅ Xử lý điều kiện biên đúng
- ✅ Tính toán phân tán (MPI_Reduce cho checksum)
- ✅ Đo thời gian bằng MPI_Wtime
- ✅ Giải phóng bộ nhớ và MPI_Finalize

## 10. Kết luận

Thiết kế song song cho LBM tương đối đơn giản nhờ:
- Tính cục bộ cao của thuật toán.
- Cấu trúc lưới đều đặn.
- Lượng giao tiếp nhỏ (chỉ 2 cột biên).

Với thiết kế hợp lý, có thể đạt được speedup gần tuyến tính khi số tiến trình không quá lớn so với kích thước lưới.
