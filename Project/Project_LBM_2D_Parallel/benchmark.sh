#!/bin/bash

echo "======================================"
echo "LBM BENCHMARK SCRIPT"
echo "======================================"

echo ""
echo "Dang bien dich..."
gcc -O2 -Wall -o lbm_serial src/lbm_serial.c -lm
mpicc -O2 -Wall -o lbm_mpi src/lbm_mpi.c -lm

if [ $? -ne 0 ]; then
    echo "❌ Loi bien dich!"
    exit 1
fi

echo "✅ Bien dich thanh cong"

CONFIG_FILE="config.txt"

if [ ! -f "$CONFIG_FILE" ]; then
    echo "nx=256" > "$CONFIG_FILE"
    echo "ny=64" >> "$CONFIG_FILE"
    echo "nsteps=10000" >> "$CONFIG_FILE"
    echo "omega=1.0" >> "$CONFIG_FILE"
    echo "u0=0.1" >> "$CONFIG_FILE"
fi

RESULTS_LOG="benchmark_results.txt"
: > "$RESULTS_LOG"

echo ""
echo "======================================"
echo "CHAY SERIAL (P=1 baseline)"
echo "======================================"
echo "--- BENCHMARK_P=1 ---" >> "$RESULTS_LOG"
./lbm_serial -c "$CONFIG_FILE" | tee -a "$RESULTS_LOG"

for np in 2 4 8; do
    echo ""
    echo "======================================"
    echo "CHAY MPI VOI $np TIEN TRINH"
    echo "======================================"
    echo "--- BENCHMARK_P=$np ---" >> "$RESULTS_LOG"
    mpirun -np $np ./lbm_mpi -c "$CONFIG_FILE" | tee -a "$RESULTS_LOG"
done

echo ""
echo "======================================"
echo "DANG XU LY KET QUA VA VE BIEU DO..."
echo "======================================"
python3 visualize_results.py
python3 plot_performance.py "$RESULTS_LOG"

echo ""
echo "======================================"
echo "HOAN TAT TOAN BO QUY TRINH"
echo "======================================"
