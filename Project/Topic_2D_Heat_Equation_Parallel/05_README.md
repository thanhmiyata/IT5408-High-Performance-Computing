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

3. Bien dich (vi du)
- Serial:
  gcc -O2 -o heat_serial src/heat_serial.c -lm
- MPI:
  mpicc -O2 -o heat_mpi src/heat_mpi.c -lm

4. Chay (vi du)
- Serial:
  ./heat_serial n K c
- MPI:
  mpirun -np 4 ./heat_mpi n K c

Trong do:
- n: so diem noi bo moi chieu (luoi n+2).
- K: so buoc thoi gian.
- c: he so khuech tan.

5. Kiem tra
- Kiem tra dieu kien on dinh: delta_t <= delta_s^2 / (2c).
- Neu can, in ra gia tri u tai mot vai diem de so sanh.

6. Do hieu nang
- Thoi gian chay lay bang MPI_Wtime (rank 0).
- Speedup = T1 / Tp, Efficiency = Speedup / p.
