#!/usr/bin/env python3

import re
import sys
import matplotlib.pyplot as plt
import numpy as np


def parse_benchmark_log(filepath):
    """Đọc thời gian thực thi và MLUPS thực tế từ file log của benchmark.sh."""
    with open(filepath, "r", encoding="utf-8") as f:
        content = f.read()

    sections = re.split(r"--- BENCHMARK_P=(\d+) ---", content)
    procs, times, mlups = [], [], []
    for idx in range(1, len(sections), 2):
        p = int(sections[idx])
        body = sections[idx + 1]
        t_match = re.search(r"Thoi gian tinh:\s*([\d.]+)\s*giay", body)
        m_match = re.search(r"MLUPS:\s*([\d.]+)", body)
        if t_match and m_match:
            procs.append(p)
            times.append(float(t_match.group(1)))
            mlups.append(float(m_match.group(1)))

    if not procs:
        raise ValueError(f"Khong tim thay du lieu benchmark hop le trong {filepath}")

    return np.array(procs), np.array(times), np.array(mlups)


def main():
    if len(sys.argv) > 1:
        log_file = sys.argv[1]
        print(f"Dang doc so lieu thuc nghiem tu: {log_file}")
        processes, times, mlups = parse_benchmark_log(log_file)
    else:
        print("Khong truyen file log -> su dung so lieu mac dinh.")
        processes = np.array([1, 2, 4, 8])
        times = np.array([1.580, 0.900, 0.539, 0.950])
        mlups = np.array([103.70, 182.04, 303.71, 172.46])

    time_serial = times[0]
    speedup = time_serial / times
    efficiency = (speedup / processes) * 100
    ideal_speedup = processes

    fig, axes = plt.subplots(1, 3, figsize=(15, 4))

    ax1 = axes[0]
    ax1.plot(processes, ideal_speedup, 'k--', label='Lý tưởng', linewidth=2)
    ax1.plot(processes, speedup, 'ro-', label='Thực tế', linewidth=2, markersize=8)
    ax1.set_xlabel('Số tiến trình (P)', fontsize=12)
    ax1.set_ylabel('Speedup', fontsize=12)
    ax1.set_title('Speedup vs Số tiến trình', fontsize=14, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    ax1.legend(fontsize=10)
    ax1.set_xticks(processes)

    for p, s in zip(processes, speedup):
        ax1.annotate(f'{s:.2f}×',
                     xy=(p, s),
                     xytext=(5, 5),
                     textcoords='offset points',
                     fontsize=9)

    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b']
    ax2 = axes[1]
    ax2.bar(processes, mlups, color=colors[:len(processes)],
            alpha=0.75, edgecolor='black', linewidth=1.5)
    ax2.set_xlabel('Số tiến trình (P)', fontsize=12)
    ax2.set_ylabel('MLUPS', fontsize=12)
    ax2.set_title('MLUPS vs Số tiến trình', fontsize=14, fontweight='bold')
    ax2.grid(True, alpha=0.3, axis='y')
    ax2.set_xticks(processes)

    for p, m in zip(processes, mlups):
        ax2.text(p, m + max(mlups) * 0.03, f'{m:.1f}',
                 ha='center', va='bottom', fontsize=10, fontweight='bold')

    ax3 = axes[2]
    ax3.plot([processes.min(), processes.max()], [100, 100], 'k--', label='Lý tưởng (100%)', linewidth=2)
    ax3.plot(processes, efficiency, 'gs-', label='Thực tế', linewidth=2, markersize=8)
    ax3.set_xlabel('Số tiến trình (P)', fontsize=12)
    ax3.set_ylabel('Efficiency (%)', fontsize=12)
    ax3.set_title('Efficiency vs Số tiến trình', fontsize=14, fontweight='bold')
    ax3.grid(True, alpha=0.3)
    ax3.legend(fontsize=10)
    ax3.set_xticks(processes)
    ax3.set_ylim([0, 115])

    for p, e in zip(processes, efficiency):
        ax3.annotate(f'{e:.1f}%',
                     xy=(p, e),
                     xytext=(5, -15),
                     textcoords='offset points',
                     fontsize=9)

    plt.tight_layout()
    plt.savefig('performance_charts.png', dpi=300, bbox_inches='tight')
    print("✅ Đã lưu biểu đồ vào: performance_charts.png")

    print("\n" + "=" * 60)
    print("BẢNG KẾT QUẢ HIỆU NĂNG")
    print("=" * 60)
    print(f"{'P':>3} | {'Time (s)':>10} | {'MLUPS':>8} | {'Speedup':>8} | {'Efficiency':>10}")
    print("-" * 60)
    for i in range(len(processes)):
        print(f"{processes[i]:>3} | {times[i]:>10.3f} | {mlups[i]:>8.2f} | {speedup[i]:>8.2f}× | {efficiency[i]:>9.1f}%")
    print("=" * 60)


if __name__ == '__main__':
    main()
