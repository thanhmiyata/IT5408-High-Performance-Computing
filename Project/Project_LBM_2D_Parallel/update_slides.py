#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script cap nhat Slides 0402.pdf:
1. Doi 'Nhom 4' thanh 'Nhom 1' tren trang tieu de (Page 1)
2. Tao slide moi ve Danh sach 6 thanh vien Nhom 1 chuan giao dien Beamer
3. Chen slide thanh vien vao ngay sau trang tieu de.
"""

import os
import subprocess
import pymupdf

def create_member_slide_pdf():
    current_dir = os.path.dirname(os.path.abspath(__file__))
    html_file = os.path.join(current_dir, "slide_members.html")
    pdf_file = os.path.join(current_dir, "slide_members.pdf")

    html_content = """<!DOCTYPE html>
<html>
<head>
<meta charset="UTF-8">
<style>
@page {
    size: 160mm 90mm; /* 453.54pt x 255.12pt (16:9 Beamer standard) */
    margin: 0;
}
body {
    margin: 0;
    padding: 0;
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif;
    background-color: #ffffff;
    color: #1a1a1a;
    width: 160mm;
    height: 90mm;
    box-sizing: border-box;
}
.header-bar {
    background-color: #0b3c5d;
    color: #ffffff;
    padding: 6mm 10mm;
    display: flex;
    justify-content: space-between;
    align-items: center;
}
.header-title {
    font-size: 13pt;
    font-weight: 700;
    letter-spacing: 0.5px;
}
.header-sub {
    font-size: 8.5pt;
    opacity: 0.9;
}
.content {
    padding: 5mm 10mm;
}
table {
    width: 100%;
    border-collapse: collapse;
    font-size: 8pt;
    margin-top: 2mm;
}
th {
    background-color: #f0f4f8;
    color: #0b3c5d;
    padding: 4px 8px;
    text-align: left;
    border-bottom: 2px solid #0b3c5d;
    font-weight: 700;
}
td {
    padding: 3.5px 8px;
    border-bottom: 1px solid #e2e8f0;
}
tr:nth-child(even) {
    background-color: #f8fafc;
}
.role {
    color: #2b6cb0;
    font-weight: 600;
}
.footer {
    position: absolute;
    bottom: 3mm;
    left: 10mm;
    right: 10mm;
    display: flex;
    justify-content: space-between;
    font-size: 6.5pt;
    color: #718096;
    border-top: 0.5px solid #e2e8f0;
    padding-top: 1.5mm;
}
</style>
</head>
<body>
<div class="header-bar">
    <div class="header-title">Thành viên Thực hiện - Nhóm 1</div>
    <div class="header-sub">Học phần IT5408: Tính toán hiệu năng cao</div>
</div>

<div class="content">
    <table>
        <thead>
            <tr>
                <th style="width: 8%; text-align: center;">STT</th>
                <th style="width: 28%;">Họ và tên</th>
                <th style="width: 18%;">Mã học viên</th>
                <th style="width: 46%;">Phân công nhiệm vụ chính</th>
            </tr>
        </thead>
        <tbody>
            <tr>
                <td style="text-align: center;"><b>1</b></td>
                <td><b>Phạm Tiến Thành</b></td>
                <td><b>20242422M</b></td>
                <td class="role">Trưởng nhóm, Cài đặt song song MPI, Benchmark</td>
            </tr>
            <tr>
                <td style="text-align: center;">2</td>
                <td>Trần Thị Hạnh</td>
                <td>20242058M</td>
                <td>Phân tích mô hình toán D2Q9, Điều kiện biên vật lý</td>
            </tr>
            <tr>
                <td style="text-align: center;">3</td>
                <td>Bùi Minh Đức</td>
                <td>20242056M</td>
                <td>Cài đặt thuật toán tuần tự C99, Kiểm thử Checksum</td>
            </tr>
            <tr>
                <td style="text-align: center;">4</td>
                <td>Nguyễn Thị Huê</td>
                <td>20252584M</td>
                <td>Xử lý dữ liệu mô phỏng, Trực quan hóa kết quả</td>
            </tr>
            <tr>
                <td style="text-align: center;">5</td>
                <td>Nguyễn Đình Hiếu</td>
                <td>20261223M</td>
                <td>Thiết kế Slide Beamer, Phân tích Scaling hiệu năng</td>
            </tr>
            <tr>
                <td style="text-align: center;">6</td>
                <td>Đỗ Minh Hiếu</td>
                <td>20242031M</td>
                <td>Biên soạn báo cáo kỹ thuật, Phân tích Overhead</td>
            </tr>
        </tbody>
    </table>
</div>

<div class="footer">
    <span>Trường Công nghệ Thông tin và Truyền thông - Đại học Bách Khoa Hà Nội</span>
    <span>GVHD: TS. Nguyễn Quốc Tuấn</span>
</div>
</body>
</html>
"""

    with open(html_file, "w", encoding="utf-8") as f:
        f.write(html_content)

    chrome_path = "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"
    cmd = [
        chrome_path,
        "--headless",
        "--disable-gpu",
        "--no-sandbox",
        "--run-all-compositor-stages-before-draw",
        "--virtual-time-budget=3000",
        f"--print-to-pdf={pdf_file}",
        "--no-pdf-header-footer",
        f"file://{html_file}"
    ]
    subprocess.run(cmd, check=True)
    return pdf_file

def update_slides_pdf():
    current_dir = os.path.dirname(os.path.abspath(__file__))
    slides_pdf = os.path.join(current_dir, "Slides 0402.pdf")

    # 1. Tao slide thanh vien PDF
    print("=> Dang tao slide Danh sach thanh vien Nhom 1...")
    member_pdf_path = create_member_slide_pdf()

    # 2. Mo doc cu
    doc = pymupdf.open(slides_pdf)
    print(f"=> Mo file {slides_pdf} (so trang goc: {len(doc)})")

    # 3. Redact 'Nhóm 4' tren trang tieu de (Page 1)
    page_title = doc[1]
    # Toa do '4' tren trang 1
    rect_4 = pymupdf.Rect(49.0, 180.0, 60.0, 193.0)
    page_title.add_redact_annot(rect_4, text="1", fill=(1, 1, 1))
    page_title.apply_redactions()
    print("=> Da cap nhat 'Nhom 4' -> 'Nhom 1' tren trang tieu de.")

    # 4. Chen slide danh sach thanh vien vao vi tri index 2 (ngay sau trang tieu de)
    member_doc = pymupdf.open(member_pdf_path)
    doc.insert_pdf(member_doc, from_page=0, to_page=0, start_at=2)
    print(f"=> Da chen slide thanh vien vao trang 2 (tong so trang moi: {len(doc)}).")

    # 5. Luu de file
    output_pdf = os.path.join(current_dir, "Slides_0402_updated.pdf")
    doc.save(output_pdf)
    doc.close()
    member_doc.close()

    # Thay the file goc
    os.replace(output_pdf, slides_pdf)
    if os.path.exists(member_pdf_path):
        os.remove(member_pdf_path)
    html_f = os.path.join(current_dir, "slide_members.html")
    if os.path.exists(html_f):
        os.remove(html_f)

    print(f"=> HOAN TAT CAP NHAT FILE SLIDE: {slides_pdf}")

if __name__ == "__main__":
    update_slides_pdf()
