# KIẾN THỨC CẦN HIỂU VỀ LATTICE BOLTZMANN METHOD

## 1. Khái niệm cơ bản

### 1.1. LBM là gì?
**Lattice Boltzmann Method (LBM)** là phương pháp số để mô phỏng dòng chảy chất lưu dựa trên lý thuyết động học phân tử (kinetic theory).

**Khác biệt chính:**
- **FDM/FEM:** Giải trực tiếp phương trình Navier-Stokes (phương trình vi phân đạo hàm riêng).
- **LBM:** Giải phương trình Boltzmann (mô tả phân bố hạt vi mô), sau đó suy ra các đại lượng vĩ mô (vận tốc, mật độ, áp suất).

### 1.2. Tại sao dùng LBM?
- **Đơn giản hóa:** Thay vì giải hệ phương trình phi tuyến phức tạp, LBM chỉ cần 2 bước: **Collision** và **Streaming**.
- **Song song tự nhiên:** Mỗi điểm lưới độc lập với nhau trong bước streaming.
- **Xử lý biên dễ dàng:** Phương pháp bounce-back cho phép mô phỏng vật cản phức tạp.

## 2. Mô hình D2Q9

### 2.1. Ý nghĩa
- **D2:** 2 chiều (2D)
- **Q9:** 9 hướng vận tốc (9 directions)

### 2.2. Các hướng vận tốc
Tại mỗi điểm lưới, có 9 hàm phân bố f₀, f₁, ..., f₈ tương ứng với 9 hướng:

```
Sơ đồ hướng (i = 0..8):
  6   2   5
    ↖ ↑ ↗
  3 ← 0 → 1
    ↙ ↓ ↘
  7   4   8

Tọa độ vận tốc:
i=0: (0, 0)   - Đứng yên
i=1: (1, 0)   - Phải
i=2: (0, 1)   - Trên
i=3: (-1, 0)  - Trái
i=4: (0, -1)  - Dưới
i=5: (1, 1)   - Phải-Trên
i=6: (-1, 1)  - Trái-Trên
i=7: (-1, -1) - Trái-Dưới
i=8: (1, -1)  - Phải-Dưới
```

### 2.3. Trọng số (Weights)
Mỗi hướng có trọng số w[i]:
- w[0] = 4/9  (hướng đứng yên)
- w[1,2,3,4] = 1/9  (4 hướng chính)
- w[5,6,7,8] = 1/36 (4 hướng chéo)

## 3. Thuật toán LBM

### 3.1. Hai bước chính

#### **Bước 1: Collision (Va chạm)**
Tính hàm phân bố sau va chạm f*ᵢ:

```
f*ᵢ(x,t) = fᵢ(x,t) - ω[fᵢ(x,t) - fᵢᵉᑫ(x,t)]
```

Trong đó:
- `ω = 1/τ`: Tham số thư giãn (relaxation parameter)
- `τ`: Thời gian thư giãn (liên quan đến độ nhớt)
- `fᵢᵉᑫ`: Hàm phân bố cân bằng (equilibrium distribution)

#### **Bước 2: Streaming (Lan truyền)**
Di chuyển hàm phân bố theo hướng vận tốc:

```
fᵢ(x + cᵢ, t+1) = f*ᵢ(x,t)
```

Trong đó `cᵢ` là vector vận tốc của hướng i.

### 3.2. Hàm phân bố cân bằng
```
fᵢᵉᑫ = wᵢ × ρ × [1 + 3(cᵢ·u) + 9/2(cᵢ·u)² - 3/2(u·u)]
```

Trong đó:
- `ρ`: Mật độ
- `u`: Vận tốc vĩ mô
- `cᵢ`: Vector vận tốc rời rạc hướng i
- `wᵢ`: Trọng số hướng i

### 3.3. Tính các đại lượng vĩ mô
Từ hàm phân bố, ta tính:

**Mật độ:**
```
ρ = Σ fᵢ  (tổng trên tất cả 9 hướng)
```

**Vận tốc:**
```
ρ × u = Σ (fᵢ × cᵢ)
```

## 4. Điều kiện biên

### 4.1. Bounce-back (Tường cứng)
Khi hạt va chạm với tường, nó bật ngược lại:
```
fᵢ(x,t+1) = f_opposite(x,t)
```

Ví dụ: Nếu hạt đi sang phải (i=1) va tường, nó sẽ bật về trái (i=3).

### 4.2. Zou-He (Biên vận tốc/áp suất)
Dùng để đặt vận tốc hoặc mật độ tại biên đầu vào/ra.

### 4.3. Periodic (Chu kỳ)
Biên trái nối với biên phải (dùng cho dòng chảy tuần hoàn).

## 5. Tham số vật lý

### 5.1. Độ nhớt động học
```
ν = (τ - 0.5) / 3
```

### 5.2. Số Reynolds
```
Re = (u₀ × L) / ν
```

Trong đó:
- `u₀`: Vận tốc đặc trưng
- `L`: Chiều dài đặc trưng
- `ν`: Độ nhớt động học

### 5.3. Điều kiện ổn định
Để mô phỏng ổn định:
- `τ > 0.5` (thường chọn τ = 0.6 - 1.0)
- `u₀ < 0.1` (vận tốc phải nhỏ để đảm bảo tính nén được thấp)

## 6. Quy trình tính toán

```
1. Khởi tạo:
   - Đặt ρ = 1.0, u = (0,0) tại mọi điểm
   - Tính fᵢ = fᵢᵉᑫ

2. Vòng lặp thời gian (t = 0 → T):
   a. Collision: Tính f*ᵢ từ fᵢ và fᵢᵉᑫ
   b. Streaming: Di chuyển f*ᵢ theo hướng cᵢ
   c. Áp dụng điều kiện biên
   d. Tính ρ và u từ fᵢ mới
   e. (Tùy chọn) Ghi kết quả ra file

3. Xuất kết quả cuối cùng
```

## 7. Ưu điểm của LBM

1. **Đơn giản:** Chỉ có phép cộng, trừ, nhân (không có đạo hàm phức tạp).
2. **Song song cao:** Mỗi điểm lưới độc lập trong bước streaming.
3. **Xử lý biên dễ:** Bounce-back rất đơn giản để code.
4. **Mở rộng tốt:** Dễ thêm lực ngoài, nhiều pha, nhiệt độ.

## 8. Nhược điểm của LBM

1. **Bộ nhớ lớn:** Cần lưu 9 giá trị fᵢ tại mỗi điểm (so với FDM chỉ cần 1-2 giá trị).
2. **Giới hạn vận tốc:** u phải nhỏ (< 0.1) để đảm bảo tính nén được thấp.
3. **Khó hiểu:** Lý thuyết động học phân tử phức tạp hơn so với phương trình vi phân.

## 9. Ứng dụng thực tế

- **Khí động học:** Thiết kế ô tô, máy bay
- **Vi lưu:** Dòng chảy trong chip, lab-on-a-chip
- **Địa chất:** Dòng chảy trong đá xốp
- **Y sinh:** Dòng máu trong mạch
- **Môi trường:** Phát tán khí thải

## 10. Tài liệu tham khảo

1. **Tutorial cơ bản:** https://www.youtube.com/watch?v=ZUXmO4hu-20
2. **Paper gốc:** S. Chen & G.D. Doolen, "Lattice Boltzmann Method for Fluid Flows" (1998)
3. **Sách:** T. Krüger et al., "The Lattice Boltzmann Method" (2017)
