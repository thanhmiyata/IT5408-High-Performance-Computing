Mô hình toán học (2D Heat Equation)

1. Bài toán
u_t = c (u_xx + u_yy), 0 <= x,y <= 1, t >= 0

Điều kiện ban đầu:
u(0, x, y) = f(x, y)

Điều kiện biên:
u(t, 0, y) = alpha0(y)
u(t, 1, y) = alpha1(y)
u(t, x, 0) = beta0(x)
u(t, x, 1) = beta1(x)

2. Rời rạc hóa
delta_s = 1/(n+1)
x_i = i * delta_s, y_j = j * delta_s, i,j = 0..n+1
delta_t: chọn theo điều kiện ổn định
t_k = k * delta_t

3. Sai phân
u_t ~ (u^{k+1}_{i,j} - u^k_{i,j}) / delta_t
u_xx ~ (u^k_{i+1,j} - 2u^k_{i,j} + u^k_{i-1,j}) / delta_s^2
u_yy ~ (u^k_{i,j+1} - 2u^k_{i,j} + u^k_{i,j-1}) / delta_s^2

4. Công thức cập nhật (explicit)
r = c * delta_t / delta_s^2
u^{k+1}_{i,j} = u^k_{i,j} + r * (u^k_{i+1,j} + u^k_{i-1,j}
                                + u^k_{i,j+1} + u^k_{i,j-1}
                                - 4*u^k_{i,j})

5. Điều kiện ổn định
delta_t <= (delta_s^2) / (2*c)

6. Đầu vào và đầu ra
- Đầu vào: n, delta_t, c, số bước thời gian K, f(x,y), các điều kiện biên.
- Đầu ra: u^K (trường nhiệt độ tại thời điểm cuối).
