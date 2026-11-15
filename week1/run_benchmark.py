#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script Benchmark thuc nghiem theo huong dan cua Giang vien:
- Do Runtime dem so nguyen to theo N (1 trieu, 10 trieu, 50 trieu, 100 trieu)
- Chay tren 1, 2, 4, 8 luong (Apple M3 - 8 cores)
- Tinh toan Speedup thuc te va Potential Speedup theo Dinh luat Amdahl
- Xuat ra bang Markdown va bieu do truc quan
"""

import subprocess
import time
import os
import matplotlib.pyplot as plt
import numpy as np

N_LIST = [
    (1000000, "1 Triệu (10^6)"),
    (10000000, "10 Triệu (10^7)"),
    (50000000, "50 Triệu (5x10^7)"),
    (100000000, "100 Triệu (10^8)")
]

THREADS = [1, 2, 4, 8]

def run_prime(n, threads):
    cmd = ["./prime_count", str(n), str(threads)]
    out = subprocess.check_output(cmd, text=True)
    # Parse "Thoi gian: X.XXXX giay."
    for line in out.splitlines():
        if "Thoi gian:" in line:
            parts = line.split()
            sec = float(parts[2])
            return sec
    return 0.0

def main():
    print("=" * 70)
    print("THUC NGHIEM DO HIỆU NĂNG THEO ĐỊNH LUẬT AMDAHL & SỐ LUỒNG (1, 2, 4, 8)")
    print("=" * 70)

    results = []

    for n, label in N_LIST:
        print(f"\n>>> Kich thuoc N = {label} <<<")
        # Chay voi 1 luong (tuan tu)
        t1 = run_prime(n, 1)
        print(f"  P = 1 luong (Serial): T1 = {t1:.4f}s")
        
        row = {
            "n": n,
            "label": label,
            "times": {1: t1},
            "speedup": {1: 1.0},
            "efficiency": {1: 1.0}
        }

        for p in [2, 4, 8]:
            tp = run_prime(n, p)
            sp = t1 / tp if tp > 0 else 0
            eff = sp / p
            row["times"][p] = tp
            row["speedup"][p] = sp
            row["efficiency"][p] = eff
            print(f"  P = {p} luong: Tp = {tp:.4f}s | Speedup = {sp:.2f}x | Efficiency = {eff*100:.1f}%")
            
        results.append(row)

    # Tinh toan Potential Speedup theo Amdahl cho N = 100M
    # Gia su ty le song song hoa thuc nghiem uoc tinh tai 4 luong:
    # sp = 1 / ((1-p_frac) + p_frac/4) => (1-p_frac) + p_frac/4 = 1/sp
    # 1 - p_frac * (3/4) = 1/sp => p_frac = (1 - 1/sp) * (4/3)
    target_row = results[-1] # N = 100M
    sp4 = target_row["speedup"][4]
    p_estimated = min(0.999, max(0.8, (1.0 - 1.0 / sp4) * (4.0 / 3.0)))
    s_estimated = 1.0 - p_estimated

    print("\n" + "=" * 70)
    print(f"PHAN TICH AMDAHL CHO N = 100 TRIEU:")
    print(f"  Uoc luong phan song song thuc te (p): {p_estimated*100:.2f}%")
    print(f"  Uoc luong phan tuan tu (s):            {s_estimated*100:.2f}%")
    print(f"  Speedup toi da ly thuyet (N -> inf):   S_max = 1/s = {1.0/s_estimated:.2f}x")
    print("=" * 70)

    # In bang Markdown
    print("\n### BANG SO LIEU THUC NGHIEM:")
    print("| Kích thước N | Số luồng P | Runtime (s) | Speedup thực tế | Hiệu suất E | Amdahl Potential Speedup |")
    print("| :--- | :---: | :---: | :---: | :---: | :---: |")
    for r in results:
        for p in THREADS:
            amdahl_sp = 1.0 / (s_estimated + p_estimated / p)
            print(f"| {r['label']} | {p} | {r['times'][p]:.4f} | {r['speedup'][p]:.2f}x | {r['efficiency'][p]*100:.1f}% | {amdahl_sp:.2f}x |")

    # Ve bieu do truc quan hoa
    os.makedirs("charts", exist_ok=True)
    plt.style.use('seaborn-v0_8-whitegrid' if 'seaborn-v0_8-whitegrid' in plt.style.available else 'default')
    
    # Bieu do 1: Runtime theo N va so luong
    fig, ax = plt.subplots(figsize=(9, 5), dpi=150)
    x_indices = np.arange(len(N_LIST))
    width = 0.18
    colors = ["#4C72B0", "#55A868", "#C44E52", "#8172B2"]

    for i, p in enumerate(THREADS):
        t_vals = [r["times"][p] for r in results]
        ax.bar(x_indices + (i - 1.5) * width, t_vals, width, label=f"P = {p} luồng", color=colors[i], edgecolor='black', linewidth=0.5)

    ax.set_xlabel("Kích thước bài toán (N)", fontsize=11, fontweight='bold')
    ax.set_ylabel("Runtime (giây) - Càng thấp càng tốt", fontsize=11, fontweight='bold')
    ax.set_title("Runtime đếm số nguyên tố theo kích thước N và số luồng", fontsize=12, fontweight='bold')
    ax.set_xticks(x_indices)
    ax.set_xticklabels([r["label"] for r in results])
    ax.legend()
    fig.tight_layout()
    fig.savefig("charts/runtime_vs_n.png")
    plt.close(fig)

    # Bieu do 2: Speedup thuc te vs Amdahl Potential Speedup tren N = 100M
    fig, ax = plt.subplots(figsize=(8, 5), dpi=150)
    ax.plot(THREADS, THREADS, 'k--', label="Ideal Linear Speedup (Lý tưởng)", alpha=0.7)
    
    actual_sp = [target_row["speedup"][p] for p in THREADS]
    amdahl_sp = [1.0 / (s_estimated + p_estimated / p) for p in THREADS]
    
    ax.plot(THREADS, amdahl_sp, 'o-', color='#C44E52', linewidth=2, label=f"Amdahl Potential (p={p_estimated*100:.1f}%)")
    ax.plot(THREADS, actual_sp, 's-', color='#4C72B0', linewidth=2, label="Thực tế đo được (Apple M3)")

    ax.set_xlabel("Số luồng xử lý (P)", fontsize=11, fontweight='bold')
    ax.set_ylabel("Hệ số tăng tốc (Speedup)", fontsize=11, fontweight='bold')
    ax.set_title("So sánh Speedup thực tế và Amdahl Potential Speedup (N = 100M)", fontsize=12, fontweight='bold')
    ax.set_xticks(THREADS)
    ax.legend()
    fig.tight_layout()
    fig.savefig("charts/amdahl_speedup_comparison.png")
    plt.close(fig)
    print("\n=> Da ve xong 2 bieu do: charts/runtime_vs_n.png va charts/amdahl_speedup_comparison.png")

if __name__ == "__main__":
    main()
