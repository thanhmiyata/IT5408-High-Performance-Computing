# TỔNG QUAN ĐỀ TÀI

## Chủ đề đề xuất
**Lattice Boltzmann Method (LBM) - Mô phỏng Dòng chảy 2D Bằng Thuật Toán Song Song**

## Lý do chọn
- **Tính song song tự nhiên:** Mỗi điểm lưới chỉ tương tác với các điểm lân cận → dễ dàng phân chia miền.
- **Ứng dụng thực tế cao:** Khí động học ô tô, thiết kế cánh máy bay, mô phỏng dòng chảy trong vi mạch.
- **Hiệu quả trên kiến trúc song song:** Cấu trúc dữ liệu đều đặn, phù hợp với MPI, OpenMP, GPU.
- **Xử lý biên phức tạp:** Dễ dàng mô phỏng các vật cản hình dạng bất kỳ bằng phương pháp bounce-back.
- **Khác biệt với FDM:** LBM giải phương trình Boltzmann (mô tả phân bố hạt) thay vì trực tiếp giải phương trình Navier-Stokes.

## Bài toán cụ thể
**Mô phỏng dòng chảy Poiseuille (Dòng chảy giữa 2 tấm phẳng song song)**

- **Miền tính toán:** Lx × Ly (ví dụ: 512 × 128 điểm lưới)
- **Điều kiện biên:**
  - Trên/Dưới: Tường cứng (bounce-back)
  - Trái: Vận tốc đầu vào u₀
  - Phải: Áp suất đầu ra (hoặc điều kiện chu kỳ)
- **Số Reynolds:** Re = 100-1000
- **Mô hình:** D2Q9 (2D, 9 hướng vận tốc)

## Mục tiêu dự án
1. Hiểu rõ lý thuyết Lattice Boltzmann Method (LBM).
2. Cài đặt mô hình D2Q9 cho dòng chảy 2D.
3. Song song hóa bằng MPI với phương pháp chia miền 1D (chia theo cột).
4. So sánh hiệu năng giữa bản tuần tự và song song.
5. Đo lường MLUPS (Million Lattice Updates Per Second), Speedup, Efficiency.
6. Viết báo cáo (~10 trang) và slide thuyết trình (~10-12 slide).

## Công nghệ sử dụng
- **Ngôn ngữ:** C
- **Thư viện:** MPI (OpenMPI hoặc MPICH)
- **Phương pháp:** Lattice Boltzmann Method (LBM) - D2Q9 Model
- **Chiến lược song song:** Domain Decomposition 1D (Column-wise)

## Cấu trúc dự án
```
Project_LBM_2D_Parallel/
├── src/
│   ├── lbm_serial.c           # Mã nguồn tuần tự
│   └── lbm_mpi.c              # Mã nguồn song song MPI
├── data/
│   ├── config_small.txt       # Bộ test nhỏ (128x32)
│   └── config_large.txt       # Bộ test lớn (512x128)
├── 00_OVERVIEW.md             # File này
├── 01_KIEN_THUC_CAN_HIEU.md  # Kiến thức về LBM
├── 02_MO_HINH_TOAN.md         # Mô hình toán học D2Q9
├── 03_THIET_KE_SONG_SONG_MPI.md  # Thiết kế song song
├── 04_BAO_CAO_SLIDE_OUTLINE.md   # Khung báo cáo và slide
├── 05_README.md               # Hướng dẫn sử dụng
├── 06_TEST_DANH_GIA.md        # Hướng dẫn test và đánh giá
├── 07_DRAFT_REPORT.md         # Báo cáo chi tiết
└── 08_DRAFT_SLIDES.md         # Slide thuyết trình
```

## So sánh LBM vs FDM

| Tiêu chí | FDM (Finite Difference) | LBM (Lattice Boltzmann) |
|----------|------------------------|-------------------------|
| **Phương trình giải** | Navier-Stokes trực tiếp | Phương trình Boltzmann |
| **Đối tượng mô phỏng** | Trường vận tốc/áp suất | Hàm phân bố hạt |
| **Xử lý biên phức tạp** | Khó | Dễ (bounce-back) |
| **Tính song song** | Tốt | Rất tốt |
| **Độ chính xác** | Bậc 2 (thông thường) | Bậc 2 (D2Q9) |
| **Ứng dụng** | Nhiệt, khuếch tán | Dòng chảy, khí động học |

## Tiến độ hiện tại
- ✅ Đã chọn đề tài: Lattice Boltzmann Method
- ✅ Đã có tài liệu tham khảo
- 🔄 Đang thực hiện: Viết code LBM tuần tự và MPI
- ⏳ Chưa làm: Chạy thử nghiệm và đo hiệu năng
- ⏳ Chưa làm: Hoàn thiện báo cáo và slide

## Tài liệu tham khảo
1. **Sách:** "The Lattice Boltzmann Method: Principles and Practice" - Krüger et al. (2017)
2. **Paper:** "Lattice Boltzmann Method for Fluid Flows" - Chen & Doolen (1998)
3. **Tutorial:** http://www.lbmethod.org/
4. **Code mẫu:** http://www.palabos.org/
