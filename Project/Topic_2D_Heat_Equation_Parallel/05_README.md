README - 2D Heat Equation (MPI)

1. Mo ta
- Giai phuong trinh nhiet 2D bang FDM explicit.
- Song song hoa bang MPI (chia mien 1D theo hang).

2. Cau truc de xuat
- src/heat_serial.c
- src/heat_mpi.c
- data/ (neu can)
- report.pdf
- slides.pdf

3. Input data (file)
- Tao file text voi dinh dang key=value, vi du:
  n=100
  K=500
  c=0.1
  left=10
  right=40
  bottom=30
  top=50
  init=0
- Da co san:
  data/case1.txt
  data/case2.txt

4. Bien dich (vi du)
- Serial:
  gcc -O2 -o heat_serial src/heat_serial.c -lm
- MPI:
  mpicc -O2 -o heat_mpi src/heat_mpi.c -lm

5. Chay (vi du)
- Serial:
  ./heat_serial data/case1.txt
- MPI:
  mpirun -np 4 ./heat_mpi data/case1.txt

Trong do:
- n: so diem noi bo moi chieu (luoi n+2)
- K: so buoc thoi gian
- c: he so khuech tan
- left/right/bottom/top: bien nhiet
- init: gia tri ban dau ben trong

6. Kiem tra
- Kiem tra dieu kien on dinh: delta_t <= delta_s^2 / (2c).
- Neu can, in ra gia tri u tai mot vai diem de so sanh.
- So sanh checksum (sum) giua serial va MPI.

7. Do hieu nang
- Thoi gian chay lay bang MPI_Wtime (rank 0).
- Speedup = T1 / Tp, Efficiency = Speedup / p.
