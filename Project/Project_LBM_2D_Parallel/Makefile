CC = gcc
MPICC = mpicc
CFLAGS = -O2 -Wall
LDFLAGS = -lm

SRC_DIR = src
SERIAL_SRC = $(SRC_DIR)/lbm_serial.c
MPI_SRC = $(SRC_DIR)/lbm_mpi.c

SERIAL_BIN = lbm_serial
MPI_BIN = lbm_mpi

NP ?= 4
CONFIG ?= config.txt
FRAME_STEP ?= 200
GRP ?= 01

.PHONY: all clean run_serial run_mpi plot demo benchmark report package help

all: $(SERIAL_BIN) $(MPI_BIN)

$(SERIAL_BIN): $(SERIAL_SRC)
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)

$(MPI_BIN): $(MPI_SRC)
	$(MPICC) $(CFLAGS) -o $@ $< $(LDFLAGS)

run_serial: $(SERIAL_BIN)
	./$(SERIAL_BIN) -c $(CONFIG)

run_mpi: $(MPI_BIN)
	mpirun -np $(NP) ./$(MPI_BIN) -c $(CONFIG)

plot:
	python3 visualize_results.py
	python3 plot_performance.py

demo: $(MPI_BIN)
	mpirun -np $(NP) ./$(MPI_BIN) -c $(CONFIG) -f $(FRAME_STEP)
	python3 visualize_results.py --demo --out frames_256x64.js
	@echo "=> Hay mo file lbm_demo_ui.html tren trinh duyet de xem dong chay 2D!"

demo_all: $(MPI_BIN)
	mpirun -np $(NP) ./$(MPI_BIN) -c data/config_medium.txt -f 200
	python3 visualize_results.py --demo --out frames_256x64.js
	mpirun -np $(NP) ./$(MPI_BIN) -c data/config_large.txt -f 500
	python3 visualize_results.py --demo --out frames_512x128.js
	mpirun -np $(NP) ./$(MPI_BIN) -c data/config_xlarge.txt -f 1000
	python3 visualize_results.py --demo --out frames_1024x256.js
	@echo "=> Da tao du lieu cho ca 3 luoi! Mo lbm_demo_ui.html de xem!"

benchmark:
	chmod +x benchmark.sh
	./benchmark.sh

report:
	python3 export_report_pdf.py

package: report
	python3 package_submission.py $(GRP)

clean:
	rm -f $(SERIAL_BIN) $(MPI_BIN) lbm_result.dat lbm_frames.txt frames_data.js frames_*.js benchmark_results.txt lbm_flow.gif *.o

help:
	@echo "Cac lenh make ho tro:"
	@echo "  make              : Bien dich ca ban Serial va MPI"
	@echo "  make run_serial   : Chay ban tuan tu (Serial)"
	@echo "  make run_mpi NP=4 : Chay ban song song MPI voi 4 process"
	@echo "  make demo NP=4    : Chay MPI xuat frames (-f 200) va tao du lieu cho lbm_demo_ui.html"
	@echo "  make plot         : Chay python visualize truong van toc va bieu do"
	@echo "  make benchmark    : Chay tu dong benchmark.sh va ve bieu do thuc te"
	@echo "  make report       : Xuat Bao_cao_HPC_LBM.pdf tu file markdown"
	@echo "  make package GRP=1: Dong goi goi nop bai .zip theo chuan mon hoc"
	@echo "  make clean        : Xoa file thuc thi va file ket qua tam"
