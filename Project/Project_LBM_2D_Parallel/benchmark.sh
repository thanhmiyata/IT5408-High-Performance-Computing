#!/bin/bash

echo "======================================"
echo "LBM BENCHMARK SCRIPT"
echo "======================================"

echo ""
echo "Dang bien dich..."
gcc -O2 -o lbm_serial src/lbm_serial.c -lm
mpicc -O2 -o lbm_mpi src/lbm_mpi.c -lm

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

echo ""
echo "======================================"
echo "CHAY SERIAL (P=1 baseline)"
echo "======================================"
./lbm_serial -c "$CONFIG_FILE"

for np in 2 4 8; do
    echo ""
    echo "======================================"
    echo "CHAY MPI VOI $np TIEN TRINH"
    echo "======================================"
    mpirun -np $np ./lbm_mpi -c "$CONFIG_FILE"
done

echo ""
echo "======================================"
echo "DANG XU LY KET QUA VA VE BIEU DO..."
echo "======================================"
python3 visualize_results.py
python3 plot_performance.py

echo ""
echo "======================================"
echo "HOAN TAT TOAN BO QUY TRINH"
echo "======================================"
