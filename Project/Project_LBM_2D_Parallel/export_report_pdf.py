#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script chuyen doi Bao_cao_HPC_LBM.md sang file HTML va PDF
chuan format hoc thuat A4 theo yeu cau cua Giang vien mon HPC.
"""

import os
import subprocess
import markdown

def export_pdf():
    current_dir = os.path.dirname(os.path.abspath(__file__))
    md_file = os.path.join(current_dir, "Bao_cao_HPC_LBM.md")
    html_file = os.path.join(current_dir, "Bao_cao_HPC_LBM.html")
    pdf_file = os.path.join(current_dir, "Bao_cao_HPC_LBM.pdf")

    if not os.path.exists(md_file):
        print(f"Loi: Khong tim thay file {md_file}")
        return

    with open(md_file, "r", encoding="utf-8") as f:
        md_text = f.read()

    # Convert markdown sang HTML với các extension chuẩn
    html_body = markdown.markdown(md_text, extensions=['tables', 'fenced_code', 'nl2br'])

    full_html = f"""<!DOCTYPE html>
<html lang="vi">
<head>
<meta charset="UTF-8">
<title>Báo cáo Dự án HPC: LBM 2D Parallel with MPI</title>

<!-- MathJax CDN ho tro cong thuc LaTeX toan hoc -->
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
    line-height: 1.6;
    color: #1a1a1a;
    font-size: 10.5pt;
    margin: 0;
    padding: 0;
}}
h1 {{
    font-size: 16pt;
    text-align: center;
    color: #0d3b66;
    margin-bottom: 8px;
    text-transform: uppercase;
    font-weight: 700;
}}
h2 {{
    font-size: 13pt;
    color: #1d3557;
    border-bottom: 1.5px solid #0d3b66;
    padding-bottom: 4px;
    margin-top: 24px;
    margin-bottom: 12px;
    page-break-after: avoid;
}}
h3 {{
    font-size: 11.5pt;
    color: #2b2d42;
    margin-top: 16px;
    margin-bottom: 8px;
    page-break-after: avoid;
}}
h4 {{
    font-size: 10.5pt;
    color: #333333;
    margin-top: 12px;
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
    margin: 14px 0;
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
    font-size: 8.5pt;
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
img {{
    max-width: 95%;
    height: auto;
    display: block;
    margin: 12px auto;
    border: 1px solid #e0e0e0;
    border-radius: 4px;
    page-break-inside: avoid;
}}
ul, ol {{
    margin-top: 4px;
    margin-bottom: 8px;
    padding-left: 22px;
}}
li {{
    margin-bottom: 4px;
}}
strong {{
    color: #0f172a;
}}
</style>
</head>
<body>
{html_body}
</body>
</html>
"""

    with open(html_file, "w", encoding="utf-8") as f:
        f.write(full_html)
    print(f"=> Da tao xong file HTML: {html_file}")

    # Su dung Google Chrome headless tren macOS de xuat PDF
    chrome_path = "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"
    if not os.path.exists(chrome_path):
        print("Loi: Khong tim thay Google Chrome de xuat PDF.")
        return

    cmd = [
        chrome_path,
        "--headless",
        "--disable-gpu",
        "--no-sandbox",
        "--run-all-compositor-stages-before-draw",
        "--virtual-time-budget=6000",
        f"--print-to-pdf={pdf_file}",
        "--no-pdf-header-footer",
        f"file://{html_file}"
    ]

    print("=> Dang render bao cao sang file PDF...")
    subprocess.run(cmd, check=True)
    if os.path.exists(pdf_file):
        size_kb = os.path.getsize(pdf_file) / 1024
        print(f"=> XUAT THANH CONG FILE PDF: {pdf_file} ({size_kb:.1f} KB)")
    else:
        print("=> Loi: Chua tao duoc file PDF.")

if __name__ == "__main__":
    export_pdf()
