#!/usr/bin/env python3

import argparse
import base64
import json
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter


def load_lbm_data(filename):
    data = np.loadtxt(filename, skiprows=1)
    x = data[:, 0]
    y = data[:, 1]
    rho = data[:, 2]
    ux = data[:, 3]
    uy = data[:, 4]

    nx = len(np.unique(x))
    ny = len(np.unique(y))

    X = x.reshape(nx, ny)
    Y = y.reshape(nx, ny)
    RHO = rho.reshape(nx, ny)
    UX = ux.reshape(nx, ny)
    UY = uy.reshape(nx, ny)

    return X, Y, RHO, UX, UY, nx, ny


def plot_velocity_field(X, Y, UX, UY, nx, ny):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    U_magnitude = np.sqrt(UX**2 + UY**2)

    levels = np.linspace(0, max(U_magnitude.max(), 1e-6), 20)
    cf = ax1.contourf(X, Y, U_magnitude, levels=levels, cmap='jet')
    ax1.set_xlabel('X', fontsize=12)
    ax1.set_ylabel('Y', fontsize=12)
    ax1.set_title('Độ lớn vận tốc |u|', fontsize=14, fontweight='bold')
    ax1.set_aspect('equal')
    cbar1 = plt.colorbar(cf, ax=ax1)
    cbar1.set_label('|u|', fontsize=11)

    skip = max(nx // 30, 1)
    ax2.quiver(X[::skip, ::skip], Y[::skip, ::skip],
               UX[::skip, ::skip], UY[::skip, ::skip],
               U_magnitude[::skip, ::skip],
               cmap='jet', scale=3, width=0.003)
    ax2.set_xlabel('X', fontsize=12)
    ax2.set_ylabel('Y', fontsize=12)
    ax2.set_title('Trường vector vận tốc', fontsize=14, fontweight='bold')
    ax2.set_aspect('equal')

    plt.tight_layout()
    plt.savefig('velocity_field.png', dpi=300, bbox_inches='tight')
    plt.close(fig)
    print("✅ Đã lưu trường vận tốc vào: velocity_field.png")


def plot_velocity_profile(UX, ny):
    mid_x = UX.shape[0] // 2
    ux_profile = UX[mid_x, :]
    y_coords = np.arange(ny)

    fig = plt.figure(figsize=(8, 6))
    plt.plot(ux_profile, y_coords, 'b-o', linewidth=2, markersize=3, label='LBM (Mô phỏng)')
    plt.xlabel('Vận tốc u_x', fontsize=12)
    plt.ylabel('Y', fontsize=12)
    plt.title('Profile vận tốc theo chiều Y', fontsize=14, fontweight='bold')
    plt.grid(True, alpha=0.3)

    y_inner = y_coords[1:-1]
    ux_inner = ux_profile[1:-1]
    y_center = (ny - 1) / 2.0
    coeffs = np.polyfit(y_inner - y_center, ux_inner, 2)
    u_fit = np.polyval(coeffs, y_coords - y_center)

    plt.plot(u_fit, y_coords, 'r--', linewidth=2, label='Fit Parabolic (Lý thuyết)')
    plt.legend(fontsize=10)

    u_max = ux_profile.max()
    plt.text(0.02, ny * 0.85, f'u_max (LBM) = {u_max:.4f}',
             fontsize=10, bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5),
             transform=plt.gca().get_yaxis_transform())

    plt.tight_layout()
    plt.savefig('velocity_profile.png', dpi=300, bbox_inches='tight')
    plt.close(fig)
    print("✅ Đã lưu profile vận tốc vào: velocity_profile.png")


def plot_density_field(X, Y, RHO):
    fig = plt.figure(figsize=(10, 5))
    rmin, rmax = RHO.min(), RHO.max()
    if np.isclose(rmin, rmax):
        rmax = rmin + 1e-6
    levels = np.linspace(rmin, rmax, 20)
    cf = plt.contourf(X, Y, RHO, levels=levels, cmap='viridis')
    plt.xlabel('X', fontsize=12)
    plt.ylabel('Y', fontsize=12)
    plt.title('Trường mật độ ρ', fontsize=14, fontweight='bold')
    plt.colorbar(cf, label='ρ')
    plt.axis('equal')
    plt.tight_layout()
    plt.savefig('density_field.png', dpi=300, bbox_inches='tight')
    plt.close(fig)
    print("✅ Đã lưu trường mật độ vào: density_field.png")


def load_frames_file(filepath="lbm_frames.txt"):
    """Đọc chuỗi frames được xuất bởi --frames/-f từ lbm_serial hoặc lbm_mpi."""
    with open(filepath, "r", encoding="utf-8") as f:
        header = f.readline().strip().split()
        nx, ny, nsteps, interval, num_ranks = map(int, header[:5])
        u0 = float(header[5])
        omega = float(header[6])

        frames = []
        cell_count = nx * ny
        while True:
            line = f.readline()
            if not line:
                break
            line = line.strip()
            if not line.startswith("FRAME "):
                continue
            step = int(line.split()[1])
            buf = np.fromfile(f, dtype=np.float32, count=cell_count * 3, sep=" ")
            if buf.size != cell_count * 3:
                break
            arr = buf.reshape(nx, ny, 3)
            frames.append({
                "step": step,
                "rho": arr[:, :, 0].copy(),
                "ux": arr[:, :, 1].copy(),
                "uy": arr[:, :, 2].copy(),
            })

    return {
        "nx": nx,
        "ny": ny,
        "nsteps": nsteps,
        "interval": interval,
        "num_ranks": num_ranks,
        "u0": u0,
        "omega": omega,
        "frames": frames,
    }


def export_demo_js(frames_file="lbm_frames.txt", output_js="frames_data.js"):
    """Chuyển lbm_frames.txt thành file frames_data.js nén Int16 Base64 cho lbm_demo_ui.html."""
    data = load_frames_file(frames_file)
    frames = data["frames"]
    if not frames:
        raise ValueError("Khong co frame nao trong file lbm_frames.txt")

    vel_scale = 20000.0   # Độ chính xác 0.00005, dải [-1.6, 1.6]
    rho_scale = 100000.0  # Độ chính xác 0.00001 quanh rho = 1.0

    encoded_frames = []
    for fr in frames:
        # Lưu theo thứ tự x outer, y inner (row-major [nx, ny])
        ux_q = np.clip(np.round(fr["ux"] * vel_scale), -32767, 32767).astype("<i2")
        uy_q = np.clip(np.round(fr["uy"] * vel_scale), -32767, 32767).astype("<i2")
        rho_q = np.clip(np.round((fr["rho"] - 1.0) * rho_scale), -32767, 32767).astype("<i2")

        encoded_frames.append({
            "step": fr["step"],
            "ux": base64.b64encode(ux_q.tobytes()).decode("ascii"),
            "uy": base64.b64encode(uy_q.tobytes()).decode("ascii"),
            "rho": base64.b64encode(rho_q.tobytes()).decode("ascii"),
            "max_u": float(np.sqrt(fr["ux"]**2 + fr["uy"]**2).max()),
            "avg_rho": float(fr["rho"].mean()),
        })

    payload = {
        "nx": data["nx"],
        "ny": data["ny"],
        "nsteps": data["nsteps"],
        "interval": data["interval"],
        "num_ranks": data["num_ranks"],
        "u0": data["u0"],
        "omega": data["omega"],
        "vel_scale": vel_scale,
        "rho_scale": rho_scale,
        "frames": encoded_frames,
    }

    with open(output_js, "w", encoding="utf-8") as f:
        f.write("window.LBM_DATASETS = window.LBM_DATASETS || {};\n")
        f.write(f"window.LBM_DATASETS['{data['nx']}x{data['ny']}'] = ")
        json.dump(payload, f, separators=(",", ":"))
        f.write(";\n")
        f.write(f"window.LBM_DEMO_DATA = window.LBM_DATASETS['{data['nx']}x{data['ny']}'];\n")

    size_mb = os.path.getsize(output_js) / (1024 * 1024)
    print(f"✅ Đã xuất dữ liệu UI ({len(encoded_frames)} frames) vào: {output_js} ({size_mb:.2f} MB)")
    return data


def export_gif(data, output_gif="lbm_flow.gif"):
    """Xuất ảnh động GIF quá trình phát triển dòng chảy."""
    frames = data["frames"]
    nx, ny = data["nx"], data["ny"]
    X, Y = np.meshgrid(np.arange(nx), np.arange(ny), indexing="ij")

    max_u_global = max(
        float(np.sqrt(fr["ux"]**2 + fr["uy"]**2).max()) for fr in frames
    )
    max_u_global = max(max_u_global, data["u0"], 1e-4)

    fig, ax = plt.subplots(figsize=(10, 3.2))
    u0_mag = np.sqrt(frames[0]["ux"]**2 + frames[0]["uy"]**2).T
    im = ax.imshow(u0_mag, origin="lower", cmap="turbo", vmin=0, vmax=max_u_global, aspect="equal")
    cbar = plt.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    cbar.set_label("|u|")
    title = ax.set_title(f"LBM D2Q9 Poiseuille Flow — Step {frames[0]['step']}/{data['nsteps']}")
    ax.set_xlabel("X")
    ax.set_ylabel("Y")

    def update(i):
        fr = frames[i]
        umag = np.sqrt(fr["ux"]**2 + fr["uy"]**2).T
        im.set_data(umag)
        title.set_text(f"LBM D2Q9 Poiseuille Flow — Step {fr['step']}/{data['nsteps']} (MPI P={data['num_ranks']})")
        return [im, title]

    anim = FuncAnimation(fig, update, frames=len(frames), interval=120, blit=False)
    anim.save(output_gif, writer=PillowWriter(fps=10))
    plt.close(fig)
    print(f"✅ Đã xuất ảnh động GIF vào: {output_gif}")


def main():
    parser = argparse.ArgumentParser(description="Trực quan hóa kết quả mô phỏng LBM 2D")
    parser.add_argument("--demo", action="store_true", help="Xuất dữ liệu JS cho UI web")
    parser.add_argument("--out", default="frames_data.js", help="Tên file JS đầu ra cho --demo (mặc định: frames_data.js)")
    parser.add_argument("--gif", action="store_true", help="Xuất thêm ảnh động lbm_flow.gif")
    args = parser.parse_args()

    filename = "lbm_result.dat"
    if not args.demo and not args.gif:
        if os.path.exists(filename):
            print(f"Đang đọc file: {filename}...")
            X, Y, RHO, UX, UY, nx, ny = load_lbm_data(filename)
            print(f"✅ Đã đọc lưới {nx}×{ny}")

            print("\nĐang vẽ biểu đồ tĩnh...")
            plot_velocity_field(X, Y, UX, UY, nx, ny)
            plot_velocity_profile(UX, ny)
            plot_density_field(X, Y, RHO)

            print("\n" + "=" * 50)
            print("THỐNG KÊ KẾT QUẢ")
            print("=" * 50)
            print(f"Kích thước lưới: {nx} × {ny}")
            print(f"Mật độ trung bình: {RHO.mean():.6f}")
            print(f"Vận tốc u_x: min={UX.min():.6f}, max={UX.max():.6f}, avg={UX.mean():.6f}")
            print(f"Vận tốc u_y: min={UY.min():.6f}, max={UY.max():.6f}, avg={UY.mean():.6f}")
            print(f"Độ lớn vận tốc max: {np.sqrt(UX**2 + UY**2).max():.6f}")
            print("=" * 50)
        else:
            print(f"❌ Lỗi: Không tìm thấy file {filename}")
            print("   Hãy chạy ./lbm_serial hoặc mpirun -np 4 ./lbm_mpi trước để tạo file kết quả.")

    if (args.demo or args.gif) or os.path.exists("lbm_frames.txt"):
        if os.path.exists("lbm_frames.txt"):
            print("\nĐang xử lý chuỗi frames cho giao diện UI...")
            data = export_demo_js("lbm_frames.txt", args.out)
            if args.gif:
                export_gif(data, "lbm_flow.gif")
        elif args.demo or args.gif:
            print("❌ Không tìm thấy lbm_frames.txt. Hãy chạy kèm tham số -f 200 (ví dụ: mpirun -np 4 ./lbm_mpi -c config.txt -f 200)")


if __name__ == '__main__':
    main()
