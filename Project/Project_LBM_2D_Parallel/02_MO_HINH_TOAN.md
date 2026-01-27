# MÔ HÌNH TOÁN HỌC - LATTICE BOLTZMANN METHOD D2Q9

## 1. Phương trình Boltzmann rời rạc

### 1.1. Phương trình chính
Phương trình Lattice Boltzmann với xấp xỉ BGK (Bhatnagar-Gross-Krook):

```
fᵢ(x + cᵢΔt, t + Δt) = fᵢ(x, t) - ω[fᵢ(x, t) - fᵢᵉᑫ(x, t)]
```

Trong đó:
- `fᵢ(x, t)`: Hàm phân bố hạt tại vị trí x, thời điểm t, hướng i
- `cᵢ`: Vector vận tốc rời rạc hướng i
- `Δt = 1`: Bước thời gian (đơn vị lưới)
- `ω = 1/τ`: Tham số thư giãn
- `τ`: Thời gian thư giãn
- `fᵢᵉᑫ`: Hàm phân bố cân bằng

### 1.2. Tách thành 2 bước

**Bước 1: Collision (Va chạm)**
```
f*ᵢ(x, t) = fᵢ(x, t) - ω[fᵢ(x, t) - fᵢᵉᑫ(x, t)]
```

**Bước 2: Streaming (Lan truyền)**
```
fᵢ(x + cᵢ, t+1) = f*ᵢ(x, t)
```

## 2. Mô hình D2Q9

### 2.1. Vector vận tốc rời rạc
9 hướng vận tốc trong không gian 2D:

```c
// Tọa độ (cx, cy)
c₀ = (0, 0)     // Đứng yên
c₁ = (1, 0)     // Phải
c₂ = (0, 1)     // Trên
c₃ = (-1, 0)    // Trái
c₄ = (0, -1)    // Dưới
c₅ = (1, 1)     // Phải-Trên
c₆ = (-1, 1)    // Trái-Trên
c₇ = (-1, -1)   // Trái-Dưới
c₈ = (1, -1)    // Phải-Dưới
```

### 2.2. Trọng số
```c
w₀ = 4/9        // Hướng đứng yên
w₁,₂,₃,₄ = 1/9  // 4 hướng chính (cardinal)
w₅,₆,₇,₈ = 1/36 // 4 hướng chéo (diagonal)
```

Tổng trọng số: Σwᵢ = 1

## 3. Hàm phân bố cân bằng

### 3.1. Công thức Maxwell-Boltzmann rời rạc
```
fᵢᵉᑫ = wᵢ × ρ × [1 + 3(cᵢ·u) + 9/2(cᵢ·u)² - 3/2(u·u)]
```

Trong đó:
- `ρ`: Mật độ vĩ mô
- `u = (uₓ, uᵧ)`: Vận tốc vĩ mô
- `cᵢ·u`: Tích vô hướng của vector vận tốc rời rạc và vận tốc vĩ mô

### 3.2. Ý nghĩa
- Khi `u = 0`: fᵢᵉᑫ = wᵢ × ρ (phân bố đều)
- Khi `u ≠ 0`: Phân bố lệch về hướng chuyển động

## 4. Các đại lượng vĩ mô

### 4.1. Mật độ
```
ρ(x, t) = Σ fᵢ(x, t)  (tổng trên tất cả 9 hướng)
        i=0..8
```

### 4.2. Vận tốc
```
ρ(x, t) × u(x, t) = Σ fᵢ(x, t) × cᵢ
                    i=0..8
```

Suy ra:
```
u(x, t) = [Σ fᵢ(x, t) × cᵢ] / ρ(x, t)
```

### 4.3. Áp suất
Trong đơn vị lưới, áp suất liên hệ với mật độ:
```
p = ρ × c²ₛ
```

Với `c²ₛ = 1/3` là tốc độ âm thanh bình phương trong đơn vị lưới.

## 5. Điều kiện biên

### 5.1. Bounce-back (Tường cứng)
Khi hạt va chạm với tường, nó bật ngược lại:
```
fᵢ(x_wall, t+1) = f*_opposite(i)(x_wall, t)
```

Ví dụ: 
- Hạt đi sang phải (i=1) va tường → bật về trái (i=3)
- Hạt đi lên (i=2) va tường → bật xuống (i=4)

### 5.2. Periodic (Chu kỳ)
Biên trái nối với biên phải:
```
fᵢ(0, y, t) = fᵢ(Lx-1, y, t)
fᵢ(Lx, y, t) = fᵢ(1, y, t)
```

### 5.3. Zou-He (Vận tốc/Áp suất)
Dùng để đặt vận tốc hoặc mật độ tại biên đầu vào/ra.
(Phức tạp hơn, không trình bày chi tiết ở đây)

## 6. Tham số vật lý

### 6.1. Độ nhớt động học
```
ν = (τ - 0.5) / 3
```

Hoặc ngược lại:
```
τ = 3ν + 0.5
```

### 6.2. Số Reynolds
```
Re = (u₀ × L) / ν
```

Trong đó:
- `u₀`: Vận tốc đặc trưng
- `L`: Chiều dài đặc trưng (ví dụ: chiều cao kênh)
- `ν`: Độ nhớt động học

### 6.3. Điều kiện ổn định

**Điều kiện 1:** Vận tốc phải nhỏ (tính nén được thấp)
```
u₀ < 0.1  (thường chọn u₀ ≈ 0.05 - 0.1)
```

**Điều kiện 2:** Thời gian thư giãn
```
τ > 0.5  (thường chọn τ = 0.6 - 1.0)
```

Nếu `τ` quá gần 0.5, mô phỏng sẽ không ổn định (numerical instability).

## 7. Đơn vị lưới (Lattice Units)

Trong LBM, tất cả tính toán đều ở đơn vị lưới:
- Bước lưới không gian: Δx = 1
- Bước thời gian: Δt = 1
- Tốc độ âm thanh: cₛ = 1/√3

Để chuyển đổi sang đơn vị vật lý, cần định nghĩa các hệ số tỉ lệ.

## 8. Liên hệ với phương trình Navier-Stokes

Khi khai triển Chapman-Enskog, phương trình Lattice Boltzmann hội tụ về phương trình Navier-Stokes không nén được:

```
∂u/∂t + (u·∇)u = -∇p/ρ + ν∇²u
∇·u = 0
```

Đây là cơ sở toán học chứng minh LBM có thể mô phỏng chính xác dòng chảy chất lưu.

## 9. Ví dụ tính toán

Cho lưới 256×64, vận tốc u₀ = 0.1, muốn Re = 100:

**Bước 1:** Chọn chiều dài đặc trưng L = 64 (chiều cao kênh)

**Bước 2:** Tính độ nhớt
```
ν = u₀ × L / Re = 0.1 × 64 / 100 = 0.064
```

**Bước 3:** Tính τ
```
τ = 3ν + 0.5 = 3 × 0.064 + 0.5 = 0.692
```

**Bước 4:** Tính ω
```
ω = 1/τ = 1/0.692 ≈ 1.445
```

## 10. Tóm tắt thuật toán

```
1. Khởi tạo: fᵢ = fᵢᵉᑫ với ρ=1, u=0

2. Lặp theo thời gian:
   a. Tính ρ, u từ fᵢ
   b. Tính fᵢᵉᑫ từ ρ, u
   c. Collision: f*ᵢ = fᵢ - ω(fᵢ - fᵢᵉᑫ)
   d. Streaming: fᵢ(x+cᵢ) = f*ᵢ(x)
   e. Áp dụng điều kiện biên

3. Xuất kết quả
```

## 11. Tài liệu tham khảo

1. S. Succi, "The Lattice Boltzmann Equation for Fluid Dynamics and Beyond" (2001)
2. T. Krüger et al., "The Lattice Boltzmann Method" (2017)
3. Q. Zou & X. He, "On pressure and velocity boundary conditions for the lattice Boltzmann BGK model" (1997)
