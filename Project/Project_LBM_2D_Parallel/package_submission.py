#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script dong goi goi nop bai cuoi ky mon HPC IT5408
theo chuan format yeu cau cua Giang vien TS. Nguyen Quoc Tuan.

Cau truc file zip nop bai:
  HPC_20261B_Group<Number>_LBM_MPI.zip
  Bao gom:
    1) Report: Bao_cao_HPC_LBM.pdf, Bao_cao_HPC_LBM.md
    2) Slide: Slides 0402.pdf, 08_DRAFT_SLIDES.md
    3) C code: src/lbm_serial.c
    4) Parallel code: src/lbm_mpi.c
    5) Readme: README.md, Makefile, benchmark.sh, config.txt, data/
       Kem cac script python va hinh anh minh hoa
"""

import os
import sys
import zipfile
import shutil

def package(group_number="01"):
    current_dir = os.path.dirname(os.path.abspath(__file__))
    pkg_name = f"HPC_20261B_Group{group_number}_LBM_MPI"
    zip_filename = os.path.join(current_dir, f"{pkg_name}.zip")

    # Danh sach cac file va thu muc hop le de dong goi
    items_to_pack = [
        # 1. Report
        ("Bao_cao_HPC_LBM.pdf", f"{pkg_name}/Report/Bao_cao_HPC_LBM.pdf"),
        ("Bao_cao_HPC_LBM.md", f"{pkg_name}/Report/Bao_cao_HPC_LBM.md"),
        ("TECHNICAL_REPORT.md", f"{pkg_name}/Report/TECHNICAL_REPORT.md"),

        # 2. Slide
        ("Slides 0402.pdf", f"{pkg_name}/Slide/Slides 0402.pdf"),
        ("08_DRAFT_SLIDES.md", f"{pkg_name}/Slide/08_DRAFT_SLIDES.md"),

        # 3. C code (Serial)
        ("src/lbm_serial.c", f"{pkg_name}/src/lbm_serial.c"),

        # 4. Parallel code (MPI)
        ("src/lbm_mpi.c", f"{pkg_name}/src/lbm_mpi.c"),

        # 5. Readme & Build configs
        ("README.md", f"{pkg_name}/README.md"),
        ("Makefile", f"{pkg_name}/Makefile"),
        ("benchmark.sh", f"{pkg_name}/benchmark.sh"),
        ("config.txt", f"{pkg_name}/config.txt"),
        ("data/config_small.txt", f"{pkg_name}/data/config_small.txt"),
        ("data/config_medium.txt", f"{pkg_name}/data/config_medium.txt"),
        ("data/config_large.txt", f"{pkg_name}/data/config_large.txt"),

        # Scripts & Charts
        ("plot_performance.py", f"{pkg_name}/plot_performance.py"),
        ("visualize_results.py", f"{pkg_name}/visualize_results.py"),
        ("performance_charts.png", f"{pkg_name}/performance_charts.png"),
        ("velocity_field.png", f"{pkg_name}/velocity_field.png"),
        ("velocity_profile.png", f"{pkg_name}/velocity_profile.png"),
        ("density_field.png", f"{pkg_name}/density_field.png"),
    ]

    print(f"=> Bat dau dong goi file zip: {zip_filename}")
    with zipfile.ZipFile(zip_filename, "w", zipfile.ZIP_DEFLATED) as zf:
        for src, arc in items_to_pack:
            src_path = os.path.join(current_dir, src)
            if os.path.exists(src_path):
                zf.write(src_path, arcname=arc)
                print(f"   [+] {src} -> {arc}")
            else:
                print(f"   [!] Canh bao: Khong tim thay file {src_path}")

    size_mb = os.path.getsize(zip_filename) / (1024 * 1024)
    print(f"\n=> DONG GOI HOAN TAT!")
    print(f"   File: {zip_filename}")
    print(f"   Dung luong: {size_mb:.2f} MB")
    print(f"   (Khong chua file binary hay file du lieu tam .dat)")

if __name__ == "__main__":
    grp = "01"
    if len(sys.argv) > 1:
        grp = sys.argv[1].replace("Group", "").replace("group", "")
    package(grp)
