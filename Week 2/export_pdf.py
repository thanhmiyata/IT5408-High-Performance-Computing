#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script chuyen doi Bao_cao_Bai_tap_Tuan_2.md sang file PDF
va dong goi thanh file ZIP nop bai theo chuan yeu cau cua Giang vien.
"""

import markdown
import subprocess
import os
import zipfile

def export_html_and_pdf():
    md_file = "Bao_cao_Bai_tap_Tuan_2.md"
    html_file = "Bao_cao_Bai_tap_Tuan_2.html"
    pdf_file = "Bao_cao_Bai_tap_Tuan_2.pdf"

    with open(md_file, "r", encoding="utf-8") as f:
        md_text = f.read()

    html_content = markdown.markdown(md_text, extensions=['tables', 'fenced_code', 'nl2br'])

    full_html = f"""<!DOCTYPE html>
<html lang="vi">
<head>
<meta charset="UTF-8">
<title>Báo cáo Bài tập Tuần 2 - IT5408</title>
<script>
window.MathJax = {{
  tex: {{
    inlineMath: [['$', '$'], ['\\\\(', '\\\\)']],
    displayMath: [['$$', '$$'], ['\\\\[', '\\\\]']]
  }},
  svg: {{
    fontCache: 'global'
  }}
}};
</script>
<script id="MathJax-script" async src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js"></script>

<style>
@page {{
    size: A4;
    margin: 18mm 18mm 18mm 18mm;
}}
body {{
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif;
    line-height: 1.55;
    color: #1a1a1a;
    font-size: 11pt;
    margin: 0;
    padding: 0;
}}
h1 {{
    font-size: 18pt;
    text-align: center;
    color: #0d3b66;
    margin-bottom: 5px;
    text-transform: uppercase;
}}
h2 {{
    font-size: 14pt;
    text-align: center;
    color: #1d3557;
    margin-top: 0;
    margin-bottom: 15px;
}}
h3 {{
    font-size: 12.5pt;
    color: #2b2d42;
    border-bottom: 1.5px solid #0d3b66;
    padding-bottom: 4px;
    margin-top: 22px;
    margin-bottom: 10px;
    page-break-after: avoid;
}}
h4 {{
    font-size: 11pt;
    color: #333333;
    margin-top: 14px;
    margin-bottom: 6px;
    page-break-after: avoid;
}}
hr {{
    border: none;
    border-top: 1px solid #dcdcdc;
    margin: 15px 0;
}}
table {{
    width: 100%;
    border-collapse: collapse;
    margin: 12px 0;
    font-size: 9.5pt;
    page-break-inside: avoid;
}}
th, td {{
    border: 1px solid #b0bec5;
    padding: 6px 10px;
    text-align: left;
}}
th {{
    background-color: #f0f4f8;
    color: #102a43;
    font-weight: 600;
}}
tr:nth-child(even) {{
    background-color: #fafbfc;
}}
code {{
    font-family: "SFMono-Regular", Consolas, "Liberation Mono", Menlo, Courier, monospace;
    font-size: 9pt;
    background-color: #f1f3f5;
    padding: 2px 4px;
    border-radius: 3px;
}}
pre {{
    background-color: #f8f9fa;
    border: 1px solid #e9ecef;
    padding: 10px;
    border-radius: 4px;
    font-size: 8.5pt;
    overflow-x: auto;
    page-break-inside: avoid;
}}
pre code {{
    background-color: transparent;
    padding: 0;
}}
ul, ol {{
    margin-top: 4px;
    margin-bottom: 8px;
    padding-left: 22px;
}}
li {{
    margin-bottom: 3px;
}}
strong {{
    color: #0f172a;
}}
</style>
</head>
<body>
{html_content}
</body>
</html>
"""

    with open(html_file, "w", encoding="utf-8") as f:
        f.write(full_html)
    print(f"=> Da tao xong file HTML: {html_file}")

    chrome_path = "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"
    abs_html = os.path.abspath(html_file)
    abs_pdf = os.path.abspath(pdf_file)

    cmd = [
        chrome_path,
        "--headless",
        "--disable-gpu",
        "--no-sandbox",
        "--run-all-compositor-stages-before-draw",
        "--virtual-time-budget=5000",
        f"--print-to-pdf={abs_pdf}",
        "--no-pdf-header-footer",
        f"file://{abs_html}"
    ]
    print("=> Dang khoi chay Chrome headless de render PDF...")
    subprocess.run(cmd, check=True)
    if os.path.exists(pdf_file):
        size_kb = os.path.getsize(pdf_file) / 1024
        print(f"=> XUAT THANH CONG FILE PDF: {pdf_file} ({size_kb:.1f} KB)")
    else:
        print("=> Chua tao duoc file PDF.")

def package_zip():
    zip_name = "Bai_Tap_Tuan_2_20242422M.zip"
    files_to_include = [
        "Bao_cao_Bai_tap_Tuan_2.pdf",
        "Bao_cao_Bai_tap_Tuan_2.md",
        "tsp_solver.cpp",
        "prime_count.cpp",
        "count_paths.cpp",
        "matmul_serial.cpp",
        "matmul_thread.cpp",
        "matmul_openmp.cpp",
        "matmul_mpi.cpp",
        "common.hpp",
        "benchmark_week2.py",
        "Makefile",
        "README.md",
        # Input files
        "tsp7.inp", "tsp12.inp", "tsp13.inp", "tsp14.inp", "tsp15.inp",
        "graph01.txt", "graph02.txt", "graph03.txt", "graph04.txt", "graph05.txt"
    ]

    with zipfile.ZipFile(zip_name, "w", zipfile.ZIP_DEFLATED) as zf:
        for f in files_to_include:
            if os.path.exists(f):
                zf.write(f, arcname=f"Bai_Tap_Tuan_2_20242422M/{f}")
        
        if os.path.exists("charts"):
            for root, _, files in os.walk("charts"):
                for file in files:
                    full_p = os.path.join(root, file)
                    arc_p = os.path.join("Bai_Tap_Tuan_2_20242422M", full_p)
                    zf.write(full_p, arcname=arc_p)

    zip_size_kb = os.path.getsize(zip_name) / 1024
    print(f"\n=> DA DONG GOI FILE ZIP NOP BAI: {zip_name} ({zip_size_kb:.1f} KB)")
    print("   (File ZIP chi gom PDF, ma nguon, file input va bieu do, KHONG chua binary/.venv, < 10MB)")

if __name__ == "__main__":
    export_html_and_pdf()
    package_zip()
