Thiet ke song song MPI (chia mien 1D)

1. Y tuong chia mien
- Chia luoi theo chieu x (chia theo hang).
- Moi process xu ly mot khoi hang lien tiep.
- Can 2 hang ghost (tren/duoi) de cap nhat stencil 5 diem.

2. Kich thuoc cuc bo
Giả sử so hang noi bo la n.
So process = P.
Chia do:
rows_per_proc = n / P, du cho phan du neu n % P != 0.

Moi process co:
- local_n = so hang noi bo ma process quan ly.
- Mang local: (local_n + 2) x (n + 2)
  (cong them 2 hang ghost).

3. Trao doi bien (ghost rows)
Truoc moi buoc thoi gian:
- Neu rank > 0: gui hang dau noi bo (1) cho rank-1,
  nhan ghost tren (0) tu rank-1.
- Neu rank < P-1: gui hang cuoi noi bo (local_n) cho rank+1,
  nhan ghost duoi (local_n+1) tu rank+1.

Co the dung MPI_Sendrecv hoac MPI_Isend/Irecv.

4. Vong lap cap nhat
for k in 1..K:
  exchange ghost rows
  for i in 1..local_n:
    for j in 1..n:
      update u_new[i][j] theo cong thuc FDM
  swap(u_old, u_new)

5. Dieu kien bien
- Hang bien toan cuc: i = 0 va i = n+1.
- Chi process dau/cuoi moi gan bien tuong ung.
- Cac cot bien (j=0, j=n+1) gan tuong ung tren tat ca process.

6. Gom ket qua
Sau khi xong, gom cac khoi ve rank 0:
MPI_Gatherv tren cac hang noi bo.

7. Do hieu nang
- Do thoi gian bang MPI_Wtime (rank 0).
- Tinh speedup va efficiency.

8. Loi khuyen
- Neu n nho, overhead giao tiep lon, speedup kem.
- Neu n lon, speedup tot hon.
