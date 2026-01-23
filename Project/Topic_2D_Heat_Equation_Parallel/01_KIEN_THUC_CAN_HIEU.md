Noi dung can hieu (giai thich chi tiet):

1. Phuong trinh nhiet 2D
- Doi tuong: nhiet do u(x, y, t) tren mien vuong [0,1]x[0,1] theo thoi gian.
- Phuong trinh: u_t = c * (u_xx + u_yy), c > 0 la he so khuech tan.
- Y nghia: nhiet lan truyen tu noi nhiet cao sang thap theo thoi gian.

2. Dieu kien bien va dieu kien ban dau
- Dieu kien ban dau: u(0, x, y) = f(x, y). Cho biet nhiet do ban dau.
- Dieu kien bien: u(t, 0, y) = alpha0(y), u(t, 1, y) = alpha1(y),
  u(t, x, 0) = beta0(x), u(t, x, 1) = beta1(x).
- Trong lap trinh, bien duoc gan co dinh va khong cap nhat theo thoi gian.

3. Roi rac hoa mien (luoi khong gian va thoi gian)
- Chia doan [0,1] thanh n+1 khoang, buoc luoi delta_s = 1/(n+1).
- Cac diem luoi: x_i = i * delta_s, y_j = j * delta_s, i,j = 0..n+1.
- Thoi gian: t_k = k * delta_t, k = 0,1,2,...

4. Sai phan huu han (Finite Difference Method - FDM)
- Dao ham theo thoi gian:
  u_t ~ (u^{k+1}_{i,j} - u^k_{i,j}) / delta_t.
- Dao ham bac hai theo khong gian:
  u_xx ~ (u^k_{i+1,j} - 2u^k_{i,j} + u^k_{i-1,j}) / delta_s^2
  u_yy ~ (u^k_{i,j+1} - 2u^k_{i,j} + u^k_{i,j-1}) / delta_s^2
- Thay vao phuong trinh nhiet ta co cong thuc cap nhat (explicit):
  u^{k+1}_{i,j} = u^k_{i,j}
    + r * (u^k_{i+1,j} + u^k_{i-1,j} + u^k_{i,j+1} + u^k_{i,j-1} - 4u^k_{i,j})
  voi r = c * delta_t / delta_s^2.

5. Dieu kien on dinh (CFL)
- FDM explicit on dinh neu delta_t <= (delta_s^2) / (2c).
- Neu vi pham, ket qua se dao dong va diverge.
- Trong code, chon delta_t theo cong thuc tren de an toan.

6. Kiem chung ket qua
- Kiem tra tinh doi xung: neu dieu kien bien doi xung thi nghiem doi xung.
- Kiem tra dao ham: u khong bi nan (no NaN).
- So sanh voi nghiem tieu chuan: bai toan voi dieu kien bien don gian co the co nghiem on dinh (steady-state).
- Kiem tra hoi tu: khi tang n (luoi tinh hon), ket qua on dinh hon.

7. Ly do co the song song hoa
- Cap nhat u^{k+1}_{i,j} chi can 4 diem lang gieng (i+1,i-1,j+1,j-1).
- Do do co the chia mien theo hang (y) hay cot (x).
- Chi can trao doi bien (ghost rows/cols) voi process ke ben.

8. MPI can biet
- MPI_Init, MPI_Comm_rank, MPI_Comm_size.
- MPI_Sendrecv hoac MPI_Isend/Irecv de trao doi bien.
- MPI_Barrier de dong bo.
- Do thoi gian bang MPI_Wtime.

9. Do hieu nang
- Time1: thoi gian chay 1 process.
- TimeP: thoi gian chay P process.
- Speedup = Time1 / TimeP.
- Efficiency = Speedup / P.

10. Dinh dang du lieu
- Mang 2D (n+2 x n+2) gom bien + noi bo.
- U_old va U_new, swap moi buoc thoi gian.
- Chi cap nhat diem noi bo (1..n).
