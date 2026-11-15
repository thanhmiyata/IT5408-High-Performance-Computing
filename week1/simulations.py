#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Kiem chung va Mo phong loi giai Bai tap Tuan 1 (IT5408 - Tinh toan hieu nang cao)
- Bai 1, 2: Dinh luat Amdahl
- Bai 6: Mo phong du bao thoi tiet
- Bai 8: Mo phong CRCW Arbitrary PRAM
- Bai 9: EREW PRAM Max Reduction
- Bai 11: Pointer Jumping (List Ranking)
- Bai 12: Parallel Merging Trade-off
"""

import math

def solve_bai1():
    print("=" * 60)
    print("BAI 1: DINH LUAT AMDAHL")
    print("=" * 60)
    p = 0.95
    s = 1.0 - p
    print(f"Phan song song hoa: p = {p * 100:.1f}%")
    print(f"Phan tuan tu:       s = {s * 100:.1f}%")
    
    # 32 cores
    N = 32
    speedup_32 = 1.0 / (s + p / N)
    print(f"\n1. Speedup tren {N} loi:")
    print(f"   S({N}) = 1 / ({s:.2f} + {p:.2f} / {N}) = {speedup_32:.6f}x")
    
    # N -> infinity
    speedup_inf = 1.0 / s
    print(f"\n2. Speedup khi so loi N -> vo han:")
    print(f"   S(inf) = 1 / s = 1 / {s:.2f} = {speedup_inf:.2f}x")

def solve_bai2():
    print("\n" + "=" * 60)
    print("BAI 2: TI LE TUAN TU TOI DA DE DAT SPEEDUP = 80 TREN 100 CORES")
    print("=" * 60)
    N = 100
    S_target = 80
    # S = 1 / (s + (1-s)/N) => s + (1-s)/N = 1 / S
    # s * (1 - 1/N) + 1/N = 1 / S
    # s * (N - 1) / N = 1/S - 1/N
    # s = (N/S - 1) / (N - 1)
    s_max = (N / S_target - 1.0) / (N - 1.0)
    p_min = 1.0 - s_max
    print(f"So bo xu ly N = {N}, Speedup mong muon S = {S_target}")
    print(f"Ti le phan tuan tu toi da cho phep: s_max = {s_max:.6f} ({s_max * 100:.4f}%)")
    print(f"Phan so chinh xac: s = 1 / {round(1.0 / s_max)}")
    print(f"Ti le song song toi thieu can dat:   p_min = {p_min:.6f} ({p_min * 100:.4f}%)")

def solve_bai6():
    print("\n" + "=" * 60)
    print("BAI 6: MO PHONG DU BAO THOI TIET")
    print("=" * 60)
    # Ban dau: 1x1x1 mile, 5e8 cells
    V_initial = 5e8
    # Kich thuoc giam con 0.5 x 0.5 x 0.5
    scale_factor = (1.0 / 0.5) ** 3  # = 8
    N_cells = V_initial * scale_factor
    print(f"1. So o khi giam kich thuoc (0.5x0.5x0.5):")
    print(f"   N_cells = 5e8 * 8 = {N_cells:.2e} o ({int(N_cells):,} o)")
    
    # Thoi gian mo phong: 1 tuan
    days = 7
    hours = days * 24
    minutes = hours * 60
    seconds = minutes * 60
    
    # Buoc thoi gian: 30s
    dt = 30
    N_steps = seconds // dt
    print(f"2. So buoc thoi gian mo phong (1 tuan, buoc 30s):")
    print(f"   Tong so giay = {seconds:,}s => N_steps = {seconds} / {dt} = {N_steps:,} buoc")
    
    # Phep tinh moi o
    ops_per_cell = 200
    total_ops = N_cells * N_steps * ops_per_cell
    print(f"3. Tong so phep tinh:")
    print(f"   Total Operations = {N_cells:.2e} * {N_steps} * {ops_per_cell} = {total_ops:.4e} phep tinh")
    
    # May thuc hien 10^9 ops/sec
    speed = 1e9
    time_sec = total_ops / speed
    time_hours = time_sec / 3600
    time_days = time_hours / 24
    time_months = time_days / 30.0
    print(f"4. Thoi gian chay tren may 10^9 phep tinh/giay:")
    print(f"   Thoi gian (giay): {time_sec:,.2f} giay")
    print(f"   Thoi gian (gio):  {time_hours:,.2f} gio")
    print(f"   Thoi gian (ngay): {time_days:,.2f} ngay ({time_months:.2f} thang)")

def solve_bai8():
    print("\n" + "=" * 60)
    print("BAI 8: MO PHONG CRCW ARBITRARY CHO PHEP AND")
    print("=" * 60)
    import random
    # Gia su mang A = [FALSE, TRUE]
    # Phép AND đúng: FALSE and TRUE = FALSE
    A = [False, True]
    print(f"Dau vao: A[1] = {A[0]} (FALSE), A[2] = {A[1]} (TRUE)")
    print(f"Ket qua dung theo toan hoc: {A[0]} AND {A[1]} = {A[0] and A[1]}")
    
    # Quy tac Arbitrary: chon ngau nhien 1 processor ghi thanh cong
    trials = 1000
    wrong_count = 0
    for _ in range(trials):
        # Processor 1 muon ghi A[0], Processor 2 muon ghi A[1]
        chosen_writer = random.choice([0, 1])
        RESULT = A[chosen_writer]
        if RESULT != False:
            wrong_count += 1
            
    print(f"Thu nghiem {trials} lan chay song song voi quy tac Arbitrary:")
    print(f"So lan RESULT bi ghi thanh TRUE (sai ket qua): {wrong_count}/{trials} (~{wrong_count/trials*100:.1f}%)")
    print("=> Ket luan: Thuat toan cho ket qua khong tin cay vi Arbitrary co the chon Processor co gia tri TRUE.")

def solve_bai11():
    print("\n" + "=" * 60)
    print("BAI 11: POINTER JUMPING (LIST RANKING)")
    print("=" * 60)
    # Danh sach: 2 -> 4 -> 1 -> 3, HEAD = 2
    # Cac nut: 1, 2, 3, 4
    nodes = [1, 2, 3, 4]
    link = {
        2: 4,
        4: 1,
        1: 3,
        3: 0  # NULL
    }
    
    rank = {i: 1 for i in nodes}
    next_node = {i: link[i] for i in nodes}
    
    step = 0
    print(f"Buoc {step} (Khoi tao):")
    print(f"  RANK: {rank}")
    print(f"  NEXT: {next_node}")
    
    while any(next_node[i] != 0 for i in nodes):
        step += 1
        new_rank = rank.copy()
        new_next = next_node.copy()
        
        for i in nodes:
            if next_node[i] != 0:
                nxt = next_node[i]
                new_rank[i] = rank[i] + rank[nxt]
                new_next[i] = next_node[nxt]
                
        rank = new_rank
        next_node = new_next
        
        print(f"\nBuoc {step} (Song song):")
        print(f"  RANK: {rank}")
        print(f"  NEXT: {next_node}")
        
    print(f"\n=> Thuat toan dung sau {step} buoc vi tat ca NEXT deu bang 0.")
    print("Y nghia RANK: Khoang cach den cuoi danh sach:")
    for i in [2, 4, 1, 3]:
        print(f"  Nut {i}: RANK = {rank[i]}")

if __name__ == "__main__":
    solve_bai1()
    solve_bai2()
    solve_bai6()
    solve_bai8()
    solve_bai11()
