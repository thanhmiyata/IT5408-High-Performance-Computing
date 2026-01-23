Mo hinh toan hoc (2D Heat Equation)

1. Bai toan
u_t = c (u_xx + u_yy), 0 <= x,y <= 1, t >= 0

Dieu kien ban dau:
u(0, x, y) = f(x, y)

Dieu kien bien:
u(t, 0, y) = alpha0(y)
u(t, 1, y) = alpha1(y)
u(t, x, 0) = beta0(x)
u(t, x, 1) = beta1(x)

2. Roi rac hoa
delta_s = 1/(n+1)
x_i = i * delta_s, y_j = j * delta_s, i,j = 0..n+1
delta_t: chon theo dieu kien on dinh
t_k = k * delta_t

3. Sai phan
u_t ~ (u^{k+1}_{i,j} - u^k_{i,j}) / delta_t
u_xx ~ (u^k_{i+1,j} - 2u^k_{i,j} + u^k_{i-1,j}) / delta_s^2
u_yy ~ (u^k_{i,j+1} - 2u^k_{i,j} + u^k_{i,j-1}) / delta_s^2

4. Cong thuc cap nhat (explicit)
r = c * delta_t / delta_s^2
u^{k+1}_{i,j} = u^k_{i,j} + r * (u^k_{i+1,j} + u^k_{i-1,j}
                               + u^k_{i,j+1} + u^k_{i,j-1}
                               - 4*u^k_{i,j})

5. Dieu kien on dinh
delta_t <= (delta_s^2) / (2*c)

6. Dau vao va dau ra
- Dau vao: n, delta_t, c, so buoc thoi gian K, f(x,y), cac dieu kien bien.
- Dau ra: u^K (truong nhiet do tai thoi diem cuoi).
