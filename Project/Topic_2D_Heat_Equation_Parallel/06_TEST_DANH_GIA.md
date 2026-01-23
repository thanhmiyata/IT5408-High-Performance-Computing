Huong dan test va danh gia

1. Muc tieu test
- So sanh dung/sai: serial vs MPI.
- Danh gia thoi gian chay, speedup, efficiency.

2. Dung/sai
- Chay serial va MPI cung input.
- So sanh gia tri:
  + u_mid (gia tri o tam).
  + sum (tong tat ca gia tri u).
- Hai gia tri tren can gan nhau (sai so nho).

3. Danh gia thoi gian
Chay nhieu lan va lay trung binh:
- Serial:
  ./heat_serial data/case2.txt
- MPI:
  mpirun -np 2 ./heat_mpi data/case2.txt
  mpirun -np 4 ./heat_mpi data/case2.txt
  mpirun -np 8 ./heat_mpi data/case2.txt

Ghi Time tu output.

4. Tinh speedup va efficiency
Speedup(p) = T1 / Tp
Efficiency(p) = Speedup(p) / p

5. Ky vong
- Luoi n lon thi speedup tot hon.
- Neu n nho thi overhead giao tiep lam hieu nang giam.
