#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script cap nhat toan dien thu muc Slides:
1. Cap nhat Slides/Slides.pdf:
   - Page 2: Redact 'Sinh vien: Krizpham' -> 'Nhom 1'
   - Chen slide danh sach 6 thanh vien Nhom 1 chuan giao dien Beamer 16:9
   - Page 18: Redact 'Sinh vien' va '[Ten giang vien]' -> 'Nhom thuc hien: Nhom 1' va 'TS. Nguyen Quoc Tuan'
2. Cap nhat main.tex trong Slides.zip voi thong tin Nhom 1 va 6 thanh vien.
3. Dong bo sang Project_LBM_2D_Parallel de dong goi goi nop bai.
"""

import os
import zipfile
import shutil
import subprocess
import pymupdf

BASE_DIR = "/Users/Krizpham/Thac si/IT5408-Tính toán hiệu năng cao"
SLIDES_DIR = os.path.join(BASE_DIR, "Slides")
PROJECT_DIR = os.path.join(BASE_DIR, "Project/Project_LBM_2D_Parallel")
FONT_PATH = "/System/Library/Fonts/Supplemental/Arial.ttf"

def create_member_slide_pdf():
    html_file = os.path.join(SLIDES_DIR, "slide_members_hust.html")
    pdf_file = os.path.join(SLIDES_DIR, "slide_members_hust.pdf")

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
    padding: 5.5mm 10mm;
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
    padding: 4mm 10mm;
}
table {
    width: 100%;
    border-collapse: collapse;
    font-size: 7.8pt;
    margin-top: 1.5mm;
}
th {
    background-color: #f0f4f8;
    color: #0b3c5d;
    padding: 4.5px 8px;
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
    slides_pdf = os.path.join(SLIDES_DIR, "Slides.pdf")
    print(f"=> Dang cap nhat file: {slides_pdf}")

    member_pdf = create_member_slide_pdf()

    doc = pymupdf.open(slides_pdf)
    print(f"   Tong so trang ban dau: {len(doc)}")

    # 1. Page 2 (Trang bia title)
    page_title = doc[2]
    rect_author = pymupdf.Rect(18.0, 209.0, 120.0, 224.0)
    page_title.add_redact_annot(rect_author, fill=(1, 1, 1))
    page_title.apply_redactions()
    page_title.insert_font(fontname="ArialVN", fontfile=FONT_PATH)
    page_title.insert_text(pymupdf.Point(19.8, 220.0), "Nhóm 1", fontname="ArialVN", fontsize=10.0, color=(0, 0, 0))
    print("   [+] Da sua 'Sinh vien: Krizpham' -> 'Nhom 1' tren trang tieu de.")

    # 2. Page 18 (Trang lien he)
    page_contact = doc[18]
    rect_contact_author = pymupdf.Rect(170.0, 32.0, 280.0, 46.0)
    page_contact.add_redact_annot(rect_contact_author, fill=(1, 1, 1))

    rect_lecturer = pymupdf.Rect(230.0, 59.0, 315.0, 74.0)
    page_contact.add_redact_annot(rect_lecturer, fill=(1, 1, 1))
    page_contact.apply_redactions()

    page_contact.insert_font(fontname="ArialVN", fontfile=FONT_PATH)
    page_contact.insert_text(pymupdf.Point(173.8, 43.0), "Nhóm thực hiện: Nhóm 1", fontname="ArialVN", fontsize=10.0, color=(0, 0, 0))
    page_contact.insert_text(pymupdf.Point(232.0, 70.0), "TS. Nguyễn Quốc Tuấn", fontname="ArialVN", fontsize=10.0, color=(0, 0, 0))
    print("   [+] Da sua 'Sinh vien' va '[Ten giang vien]' tren trang Lien he.")

    # 3. Chen slide danh sach thanh vien vao vi tri index 3 (ngay sau trang tieu de)
    mem_doc = pymupdf.open(member_pdf)
    doc.insert_pdf(mem_doc, from_page=0, to_page=0, start_at=3)
    print(f"   [+] Da chen slide danh sach thanh vien vao sau trang bia (Tong so trang moi: {len(doc)}).")

    # Luu file
    tmp_out = os.path.join(SLIDES_DIR, "Slides_tmp.pdf")
    doc.save(tmp_out)
    doc.close()
    mem_doc.close()
    os.replace(tmp_out, slides_pdf)

    # Don dep file tam
    if os.path.exists(member_pdf):
        os.remove(member_pdf)
    html_f = os.path.join(SLIDES_DIR, "slide_members_hust.html")
    if os.path.exists(html_f):
        os.remove(html_f)

    print(f"=> Cap nhat thanh cong Slides/Slides.pdf!")

def update_slides_zip():
    slides_zip = os.path.join(SLIDES_DIR, "Slides.zip")
    extract_dir = os.path.join(SLIDES_DIR, "latex_source")
    if os.path.exists(extract_dir):
        shutil.rmtree(extract_dir)

    print(f"=> Dang giai nen {slides_zip}...")
    with zipfile.ZipFile(slides_zip, 'r') as zf:
        zf.extractall(extract_dir)

    main_tex = os.path.join(extract_dir, "main.tex")
    with open(main_tex, "r", encoding="utf-8") as f:
        tex = f.read()

    # Thay the thong tin tac gia
    tex = tex.replace("\\author{Sinh viên: Krizpham}", "\\author{Nhóm 1}")
    tex = tex.replace("\\textbf{Sinh viên:} Krizpham", "\\textbf{Nhóm thực hiện:} Nhóm 1")
    tex = tex.replace("\\textbf{Giảng viên:} [Tên giảng viên]", "\\textbf{Giảng viên hướng dẫn:} TS. Nguyễn Quốc Tuấn")

    # Chen slide thanh vien ngay sau \husttitlepage
    member_frame = """
	% Trang bìa
	\\husttitlepage
	
	% Slide thành viên Nhóm 1
	\\begin{frame}{Thành viên Thực hiện - Nhóm 1}
		\\begin{table}
			\\centering
			\\scriptsize
			\\begin{tabular}{|c|l|c|p{6.8cm}|}
				\\hline
				\\textbf{STT} & \\textbf{Họ và tên} & \\textbf{Mã HV} & \\textbf{Phân công nhiệm vụ chính} \\\\
				\\hline
				1 & \\textbf{Phạm Tiến Thành} & \\textbf{20242422M} & Trưởng nhóm, Cài đặt song song MPI, Benchmark \\\\
				\\hline
				2 & Trần Thị Hạnh & 20242058M & Phân tích mô hình toán D2Q9, Điều kiện biên vật lý \\\\
				\\hline
				3 & Bùi Minh Đức & 20242056M & Cài đặt thuật toán tuần tự C99, Kiểm thử Checksum \\\\
				\\hline
				4 & Nguyễn Thị Huê & 20252584M & Xử lý dữ liệu mô phỏng, Trực quan hóa kết quả \\\\
				\\hline
				5 & Nguyễn Đình Hiếu & 20261223M & Thiết kế Slide Beamer, Phân tích Scaling hiệu năng \\\\
				\\hline
				6 & Đỗ Minh Hiếu & 20242031M & Biên soạn báo cáo kỹ thuật, Phân tích Overhead \\\\
				\\hline
			\\end{tabular}
		\\end{table}
	\\end{frame}
"""
    if "\\begin{frame}{Thành viên Thực hiện - Nhóm 1}" not in tex:
        tex = tex.replace("% Trang bìa\n\t\\husttitlepage", member_frame)

    with open(main_tex, "w", encoding="utf-8") as f:
        f.write(tex)

    print("   [+] Da cap nhat main.tex voi thong tin Nhom 1 va 6 thanh vien.")

    with zipfile.ZipFile(slides_zip, 'w', zipfile.ZIP_DEFLATED) as zf:
        for root, _, files in os.walk(extract_dir):
            for file in files:
                full_p = os.path.join(root, file)
                rel_p = os.path.relpath(full_p, extract_dir)
                zf.write(full_p, arcname=rel_p)

    print(f"   [+] Da dong goi lai {slides_zip}.")

    # Dong bo sang Project
    dest_beamer = os.path.join(PROJECT_DIR, "slide_beamer")
    if os.path.exists(dest_beamer):
        shutil.rmtree(dest_beamer)
    shutil.copytree(extract_dir, dest_beamer)
    print(f"   [+] Da dong bo thu muc nguon LaTeX sang: {dest_beamer}")

    shutil.copyfile(os.path.join(SLIDES_DIR, "Slides.pdf"), os.path.join(PROJECT_DIR, "Slides.pdf"))
    shutil.copyfile(os.path.join(SLIDES_DIR, "Slides.pdf"), os.path.join(PROJECT_DIR, "Slides 0402.pdf"))
    print(f"   [+] Da dong bo Slides.pdf sang Project_LBM_2D_Parallel.")

if __name__ == "__main__":
    update_slides_pdf()
    update_slides_zip()
