#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script Benchmark Toan Dien Tuan 2 (IT5408 - Tinh toan hieu nang cao):
1. Bai 1: TSP (Traveling Salesman Problem - Branch & Bound) tren tsp7, tsp12, tsp13, tsp14, tsp15
2. Bai 2: Dem so nguyen to (N = 10^6, 10^7, 5x10^7, 10^8) va doi chieu Amdahl Potential Speedup
3. Bai 3: Dem so duong di do dai k tren do thi (graph01 -> graph05)
4. Bai 4: Nhan 2 ma tran (Serial, Thread, OpenMP, MPI)
"""

import subprocess
import os
import json
import time
import matplotlib.pyplot as plt
import numpy as np

THREADS = [1, 2, 4, 8]

def run_tsp():
    print("=" * 60)
    print("1. BENCHMARK BAI TOAN 1: TSP (BRANCH & BOUND)")
    print("=" * 60)
    files = ["tsp7.inp", "tsp12.inp", "tsp13.inp", "tsp14.inp", "tsp15.inp"]
    tsp_results = []

    for f in files:
        if not os.path.exists(f): continue
        row = {"file": f, "times": {}, "speedup": {}, "min_cost": 0, "tour": ""}
        for t in THREADS:
            out = subprocess.check_output(["./tsp_solver", f, str(t)], text=True)
            for line in out.splitlines():
                if "Chi phi nho nhat" in line:
                    row["min_cost"] = int(line.split(":")[1].strip())
                elif "Chu trinh toi uu:" in line:
                    row["tour"] = line.split(":")[1].strip()
                elif "Thoi gian thuc thi:" in line:
                    parts = line.split()
                    sec = float(parts[4])
                    row["times"][t] = sec
            
            t1 = row["times"].get(1, 0.0001)
            tt = row["times"][t]
            sp = t1 / tt if tt > 0 else 1.0
            row["speedup"][t] = sp

        tsp_results.append(row)
        print(f"File {f}: Min Cost = {row['min_cost']}, Tour: {row['tour']}")
        print(f"  Times: {row['times']}")
    return tsp_results

def run_primes():
    print("\n" + "=" * 60)
    print("2. BENCHMARK BAI TOAN 2: DEM SO NGUYEN TO & DOI CHIEU AMDAHL")
    print("=" * 60)
    N_list = [
        (1000000, "1 Triệu (10^6)"),
        (10000000, "10 Triệu (10^7)"),
        (50000000, "50 Triệu (5x10^7)"),
        (100000000, "100 Triệu (10^8)")
    ]
    prime_results = []
    for n, label in N_list:
        row = {"n": n, "label": label, "times": {}, "speedup": {}, "efficiency": {}, "count": 0}
        for t in THREADS:
            out = subprocess.check_output(["./prime_count", str(n), str(t), "--json"], text=True)
            data = json.loads(out.strip())
            row["times"][t] = data["time_sec"]
            row["count"] = data["count"]
            
        t1 = row["times"][1]
        for t in THREADS:
            sp = t1 / row["times"][t]
            row["speedup"][t] = sp
            row["efficiency"][t] = sp / t
        prime_results.append(row)
        print(f"N = {label} (pi(N) = {row['count']}):")
        for t in THREADS:
            print(f"  {t} luong: {row['times'][t]:.4f}s | Speedup: {row['speedup'][t]:.2f}x | Efficiency: {row['efficiency'][t]*100:.1f}%")
    return prime_results

def run_graphs():
    print("\n" + "=" * 60)
    print("3. BENCHMARK BAI TOAN 3: DEM SO DUONG DI DO DAI K TRONG DO THI")
    print("=" * 60)
    graphs = [
        ("graph01.txt", 3),
        ("graph02.txt", 3),
        ("graph03.txt", 4),
        ("graph04.txt", 4),
        ("graph05.txt", 5)
    ]
    graph_results = []
    for f, k in graphs:
        if not os.path.exists(f): continue
        row = {"file": f, "k": k, "times": {}, "speedup": {}, "paths": 0}
        for t in THREADS:
            out = subprocess.check_output(["./count_paths", f, str(k), str(t), "--json"], text=True)
            data = json.loads(out.strip())
            row["times"][t] = data["time_sec"]
            row["paths"] = data["paths"]
            row["V"] = data["V"]
            row["E"] = data["E"]
        t1 = row["times"][1]
        for t in THREADS:
            row["speedup"][t] = t1 / row["times"][t] if row["times"][t] > 0 else 1.0
        graph_results.append(row)
        print(f"{f} (V={row['V']}, E={row['E']}, k={k}): Tong duong di = {row['paths']}")
        print(f"  Times: {row['times']}")
    return graph_results

def run_matmul():
    print("\n" + "=" * 60)
    print("4. BENCHMARK BAI TOAN BO SUNG: NHAN 2 MA TRAN (SERIAL, THREAD, OPENMP, MPI)")
    print("=" * 60)
    sizes = [(512, 512, 512), (1024, 1024, 1024)]
    matmul_results = []
    for m, n, p in sizes:
        label = f"{m}x{n}x{p}"
        # Serial
        s_out = json.loads(subprocess.check_output(["./matmul_serial", str(m), str(n), str(p), "--json"], text=True).strip())
        t_serial = s_out["time_sec"]
        
        row = {"size": label, "serial_time": t_serial, "methods": {}}
        for t in [2, 4, 8]:
            th_out = json.loads(subprocess.check_output(["./matmul_thread", str(m), str(n), str(p), str(t), "--json"], text=True).strip())
            omp_out = json.loads(subprocess.check_output(["./matmul_openmp", str(m), str(n), str(p), str(t), "--json"], text=True).strip())
            mpi_out = json.loads(subprocess.check_output(["mpirun", "-np", str(t), "./matmul_mpi", str(m), str(n), str(p), "--json"], text=True).strip())
            
            row["methods"][f"thread_{t}"] = {"time": th_out["time_sec"], "speedup": t_serial / th_out["time_sec"], "gflops": th_out["gflops"]}
            row["methods"][f"omp_{t}"] = {"time": omp_out["time_sec"], "speedup": t_serial / omp_out["time_sec"], "gflops": omp_out["gflops"]}
            row["methods"][f"mpi_{t}"] = {"time": mpi_out["time_sec"], "speedup": t_serial / mpi_out["time_sec"], "gflops": mpi_out["gflops"]}
        
        matmul_results.append(row)
        print(f"Ma tran {label} (Serial = {t_serial:.4f}s):")
        print(f"  OpenMP 4T: {row['methods']['omp_4']['time']:.4f}s | Speedup: {row['methods']['omp_4']['speedup']:.2f}x | {row['methods']['omp_4']['gflops']:.2f} GFLOPS")
        print(f"  Thread 4T: {row['methods']['thread_4']['time']:.4f}s | Speedup: {row['methods']['thread_4']['speedup']:.2f}x")
        print(f"  MPI 4P:    {row['methods']['mpi_4']['time']:.4f}s | Speedup: {row['methods']['mpi_4']['speedup']:.2f}x")
    return matmul_results

def plot_charts(tsp_res, prime_res, graph_res, matmul_res):
    os.makedirs("charts", exist_ok=True)
    plt.style.use('seaborn-v0_8-whitegrid' if 'seaborn-v0_8-whitegrid' in plt.style.available else 'default')
    
    # 1. Bieu do Prime Count: Runtime vs N & Amdahl Potential Speedup
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5), dpi=150)
    
    x = np.arange(len(prime_res))
    width = 0.2
    colors = ["#4C72B0", "#55A868", "#C44E52", "#8172B2"]
    for i, t in enumerate(THREADS):
        t_vals = [r["times"][t] for r in prime_res]
        ax1.bar(x + (i - 1.5) * width, t_vals, width, label=f"{t} luồng", color=colors[i], edgecolor='black', linewidth=0.5)
    ax1.set_xlabel("Kích thước N", fontsize=11, fontweight='bold')
    ax1.set_ylabel("Runtime (giây)", fontsize=11, fontweight='bold')
    ax1.set_title("Runtime đếm số nguyên tố theo N và số luồng", fontsize=12, fontweight='bold')
    ax1.set_xticks(x)
    ax1.set_xticklabels([r["label"] for r in prime_res], rotation=15, ha='right')
    ax1.legend()

    # Amdahl vs Actual tai N = 100M
    p_100 = prime_res[-1]
    sp4 = p_100["speedup"][4]
    p_frac = (1.0 - 1.0 / sp4) * (4.0 / 3.0)
    s_frac = 1.0 - p_frac
    amdahl_curve = [1.0 / (s_frac + p_frac / t) for t in THREADS]
    actual_curve = [p_100["speedup"][t] for t in THREADS]

    ax2.plot(THREADS, THREADS, 'k--', label="Lý tưởng (Linear)", alpha=0.7)
    ax2.plot(THREADS, amdahl_curve, 'o-', color='#C44E52', linewidth=2, label=f"Amdahl Potential (p={p_frac*100:.1f}%)")
    ax2.plot(THREADS, actual_curve, 's-', color='#4C72B0', linewidth=2, label="Thực tế đo (Apple M3)")
    ax2.set_xlabel("Số luồng xử lý (P)", fontsize=11, fontweight='bold')
    ax2.set_ylabel("Hệ số tăng tốc (Speedup)", fontsize=11, fontweight='bold')
    ax2.set_title(f"Amdahl Potential Speedup vs Thực tế (N = 100M)", fontsize=12, fontweight='bold')
    ax2.set_xticks(THREADS)
    ax2.legend()

    fig.tight_layout()
    fig.savefig("charts/prime_runtime_amdahl.png")
    plt.close(fig)

    # 2. Bieu do TSP Speedup
    fig, ax = plt.subplots(figsize=(8, 5), dpi=150)
    for r in tsp_res:
        if r["file"] in ["tsp13.inp", "tsp14.inp", "tsp15.inp"]:
            sp_vals = [r["speedup"][t] for t in THREADS]
            ax.plot(THREADS, sp_vals, marker='o', linewidth=2, label=f"{r['file']} (Cost={r['min_cost']})")
    ax.plot(THREADS, THREADS, 'k--', label="Lý tưởng (Linear)", alpha=0.6)
    ax.set_xlabel("Số luồng OpenMP (P)", fontsize=11, fontweight='bold')
    ax.set_ylabel("Hệ số tăng tốc Speedup", fontsize=11, fontweight='bold')
    ax.set_title("Hệ số tăng tốc thuật toán TSP Branch & Bound theo số luồng", fontsize=12, fontweight='bold')
    ax.set_xticks(THREADS)
    ax.legend()
    fig.tight_layout()
    fig.savefig("charts/tsp_scaling.png")
    plt.close(fig)

    # 3. Bieu do Nhan Ma tran 1024x1024x1024
    fig, ax = plt.subplots(figsize=(9, 5), dpi=150)
    m1024 = matmul_res[1]
    threads_eval = [2, 4, 8]
    x_m = np.arange(len(threads_eval))
    w = 0.25

    th_sp = [m1024["methods"][f"thread_{t}"]["speedup"] for t in threads_eval]
    omp_sp = [m1024["methods"][f"omp_{t}"]["speedup"] for t in threads_eval]
    mpi_sp = [m1024["methods"][f"mpi_{t}"]["speedup"] for t in threads_eval]

    ax.bar(x_m - w, th_sp, w, label="std::thread", color="#55A868", edgecolor='black', linewidth=0.5)
    ax.bar(x_m, omp_sp, w, label="OpenMP", color="#C44E52", edgecolor='black', linewidth=0.5)
    ax.bar(x_m + w, mpi_sp, w, label="MPI", color="#8172B2", edgecolor='black', linewidth=0.5)

    ax.set_xlabel("Số luồng / Số tiến trình", fontsize=11, fontweight='bold')
    ax.set_ylabel("Hệ số tăng tốc (Speedup so với Serial)", fontsize=11, fontweight='bold')
    ax.set_title("So sánh Speedup nhân ma trận 1024x1024x1024 (Apple M3)", fontsize=12, fontweight='bold')
    ax.set_xticks(x_m)
    ax.set_xticklabels([f"{t} Cores" for t in threads_eval])
    ax.legend()
    fig.tight_layout()
    fig.savefig("charts/matmul_comparison.png")
    plt.close(fig)
    print("=> Da ve xong toan bo bieu do vao charts/!")

def main():
    tsp_res = run_tsp()
    prime_res = run_primes()
    graph_res = run_graphs()
    matmul_res = run_matmul()
    plot_charts(tsp_res, prime_res, graph_res, matmul_res)

if __name__ == "__main__":
    main()
