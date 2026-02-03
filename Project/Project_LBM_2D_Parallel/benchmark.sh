#!/bin/bash
# Script benchmark tự động cho LBM

echo "======================================"
echo "LBM BENCHMARK SCRIPT"
echo "======================================"

# Biên dịch
echo ""
echo "Đang biên dịch..."
gcc -O2 -o lbm_serial src/lbm_serial.c -lm
mpicc -O2 -o lbm_mpi src/lbm_mpi.c -lm

if [ $? -ne 0 ]; then
    echo "❌ Lỗi biên dịch!"
    exit 1
fi

echo "✅ Biên dịch thành công"

# Tham số
NX=${1:-256}
NY=${2:-64}
NSTEPS=${3:-10000}

# Chạy serial
echo ""
echo "======================================"
echo "CHẠY SERIAL (Luoi: $NX x $NY, $NSTEPS buoc)"
echo "======================================"
./lbm_serial $NX $NY $NSTEPS

# Chạy MPI với nhiều số tiến trình
for np in 2 4 8; do
    echo ""
    echo "======================================"
    echo "CHẠY MPI VỚI $np TIẾN TRÌNH"
    echo "======================================"
    mpirun -np $np ./lbm_mpi $NX $NY $NSTEPS
done

echo ""
echo "======================================"
echo "HOÀN TẤT BENCHMARK"
echo "======================================"
echo ""
echo "Để vẽ biểu đồ hiệu năng, chạy:"
echo "  python3 plot_performance.py"
echo ""
echo "Để visualize kết quả, chạy:"
echo "  python3 visualize_results.py"
