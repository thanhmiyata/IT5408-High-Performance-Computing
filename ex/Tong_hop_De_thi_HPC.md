# Tổng hợp Đề thi Cuối kỳ IT5408 - Tính toán hiệu năng cao (HPC)

**Môn học:** IT5408 - High Performance Computing  
**Thời gian:** 60 phút  
**Định dạng:** Trắc nghiệm (Tiếng Anh) & Tự luận điền code (CUDA/C)

---

## I. PHẦN TRẮC NGHIỆM (Multiple Choice Questions)

Dưới đây là danh sách các câu hỏi và đáp án dựa trên nội dung tổng hợp từ các file ảnh:

| STT | Câu hỏi (Tiếng Anh) | Giải thích / Đáp án (Gợi ý) |
|:---:|:---|:---|
| **Q1** | Which of the following is not type of Flynn's Classification? | **B. SIDD** (Phân loại Flynn gồm: SISD, SIMD, MISD, MIMD). |
| **Q2** | In shared Memory... | **B.** Thay đổi tại một vị trí bộ nhớ bởi 1 tiến trình sẽ hiển thị với các tiến trình khác. |
| **Q3** | In distributed Memory... | **C.** Các tác vụ song song thường trao đổi dữ liệu qua mạng. |
| **Q4** | Synchronous communication operations referred to? | **B.** Phụ thuộc giữa các câu lệnh khi thứ tự thực thi ảnh hưởng đến kết quả. |
| **Q5** | Granularity is? | **A.** Tỷ lệ giữa tính toán (computation) và liên lạc (communication). |
| **Q6** | Parallel Overhead is? | **B.** Thời gian điều phối (start-up, synchronization, communication). |
| **Q7** | Scalability refers to? | **B.** Khả năng tăng tốc độ song song khi thêm nhiều bộ xử lý. |
| **Q8** | Non-Uniform Memory Access (NUMA) is? | **C.** SMP có thể truy cập bộ nhớ của bộ khác nhưng tốc độ không đồng đều. |
| **Q9** | Task dependency graph is? | **A. directed** (Đồ thị có hướng). |
| **Q10** | In the threads model... | **A.** Một tiến trình có thể chứa nhiều luồng thực thi đồng thời. |
| **Q11** | Data dependence is? | **B.** Thứ tự thực thi câu lệnh ảnh hưởng trực tiếp đến kết quả cuối cùng. |
| **Q12** | Which of followings is an example of domain decomposition? | **A. Matrix multiplication** (Phân rã miền dữ liệu cho phép nhân ma trận). |
| **Q13** | Breaking problem into discreet chunks is known as? | **A. Decomposition** (Sự phân rã bài toán). |
| **Q14** | Functional Decomposition is? | **B.** Tập trung vào các hàm/tác vụ tính toán hơn là dữ liệu. |
| **Q15** | Speedup of a program achieving 75% efficiency on 32 processors? | **B. 24** (Speedup = Efficiency * P = 0.75 * 32). |
| **Q16** | The type of parallelism supported by CUDA is best described as? | **D. MIMD** (Theo ghi chú trong ảnh đề thi). |
| **Q17** | Within a parallel region, declared variables are by default... | **B. Local** (Các biến khai báo trong vùng song song mặc định là cục bộ). |
| **Q18** | Collective communication... | **A.** Chia sẻ dữ liệu giữa tập hợp nhiều hơn 2 thành viên. |
| **Q19** | Object holding info about received message (count, tag...)? | **D. status** (Thường là `MPI_Status` trong thư viện MPI). |
| **Q20** | In a CUDA program, calling a kernel is typically referred to as? | **B. Thread invocation**. |
| **Q21** | Consider: `omp_set_num_threads(10); #pragma omp parallel { for(i=0; i<100; i++) printf(...); }`. Total prints? | **A. 1000** (10 luồng, mỗi luồng chạy vòng lặp 100 lần). |
| **Q22** | OMP_NUM_THREADS vs omp_set_num_threads overlap? | **D. omp_set_num_threads** sẽ ghi đè biến môi trường. |
| **Q23** | Which of describing correctly about GPU kernel? | **C.** Tất cả các thread block trong cùng một phép toán dùng chung một kernel. |

---

## II. PHẦN TỰ LUẬN - ĐIỀN CODE (Fill in the blanks)

**Đề bài:** Hoàn thiện chương trình CUDA tính đạo hàm bậc nhất của hàm $f(x) = \frac{8x}{x^2+1}$ bằng phương pháp Central Difference.

### 1. Central Difference Kernel (Code A)
Phần code nằm bên trong kernel `CentralDifference` để xác định vị trí của luồng:
```cuda
// --- Code A ---
int i = blockIdx.x * blockDim.x + threadIdx.x;
if (i < N) {
    // Logic tính toán đã có trong đề:
    // if (i == 0 || i == N-1) *(DF+i) = 0;
    // else *(DF+i) = (*(F+i+1) - *(F+i-1)) / (2*h);
}
```

### 2. Copy dữ liệu từ CPU sang GPU (Code B)
```cuda
// --- Code B ---
cudaMemcpy(FXgpu, FXcpu, N * sizeof(float), cudaMemcpyHostToDevice);
```

### 3. Gọi Kernel (Code C)
Cấu hình Grid và Block để thực thi trên GPU:
```cuda
// --- Code C ---
CentralDifference<<<dimGrid, dimBlock>>>(DFXgpu, FXgpu);
```

### 4. Copy kết quả từ GPU về CPU (Code D)
```cuda
// --- Code D ---
cudaMemcpy(DFXcpu, DFXgpu, N * sizeof(float), cudaMemcpyDeviceToHost);
```

---
*Tài liệu được tổng hợp tự động để hỗ trợ ôn tập cuối kỳ.*

## III. KIẾN THỨC TRỌNG TÂM CẦN HỌC THUỘC (TỦ)

Dựa trên phân tích các bài tập trong folder `ex`, dưới đây là các đoạn code "kinh điển" thường bị để trống trong đề thi:

### 1. Hàm chia việc (Workshare) cho MPI
Thường dùng trong bài Nhân ma trận để chia cột/hàng cho các process.
```c
void workshare(int rank, int nprocs, int N, int *ilo, int *ihi) {
    int nrem = N % nprocs;       // Phần dư
    int nchunk = N / nprocs;     // Kích thước chuẩn
    if (rank < nrem) {
        *ilo = 1 + rank * (nchunk + 1);
        *ihi = *ilo + nchunk;
    } else {
        *ilo = 1 + rank * nchunk + nrem;
        *ihi = *ilo + nchunk - 1;
    }
}
```

### 2. Trao đổi biên trong MPI (Halo Exchange)
Dùng trong bài toán Phương trình nhiệt 1D. Chú ý thứ tự gửi/nhận để tránh Deadlock.
```c
// Gửi biên phải cho rank+1, nhận biên trái từ rank-1
if (rank < size - 1) 
    MPI_Send(Ts + mc - 1, 1, MPI_FLOAT, rank + 1, tag, MPI_COMM_WORLD);
if (rank > 0) 
    MPI_Recv(&Tl, 1, MPI_FLOAT, rank - 1, tag, MPI_COMM_WORLD, &stat);

// Gửi biên trái cho rank-1, nhận biên phải từ rank+1
if (rank > 0) 
    MPI_Send(Ts, 1, MPI_FLOAT, rank - 1, tag, MPI_COMM_WORLD);
if (rank < size - 1) 
    MPI_Recv(&Tr, 1, MPI_FLOAT, rank + 1, tag, MPI_COMM_WORLD, &stat);
```

### 3. Chia miền tính toán trong OpenMP
Dùng khi muốn chia vòng lặp thủ công bằng `omp_get_thread_num()`.
```c
#pragma omp parallel
{
    int ID = omp_get_thread_num();
    int NT = omp_get_num_threads();
    int Mc = M / NT;
    int start = ID * Mc;
    int stop = (ID == NT - 1) ? M : (start + Mc);
    
    // Thực hiện tính toán trong khoảng [start, stop)
    #pragma omp barrier // Đợi tất cả thread tính xong biên trước khi update
}
```

### 4. Công thức tính Index trong CUDA Kernel
Bước đầu tiên và quan trọng nhất bên trong mọi Kernel.
```c
// 1D Grid/Block
int i = blockIdx.x * blockDim.x + threadIdx.x;

// 2D Grid/Block (Nếu đề bài là nhân ma trận)
int row = blockIdx.y * blockDim.y + threadIdx.y;
int col = blockIdx.x * blockDim.x + threadIdx.x;
```

### 5. Thứ tự gọi hàm MPI cơ bản
```c
MPI_Init(&argc, &argv);
MPI_Comm_rank(MPI_COMM_WORLD, &rank);
MPI_Comm_size(MPI_COMM_WORLD, &size);
// ... Tính toán ...
MPI_Finalize();
```
