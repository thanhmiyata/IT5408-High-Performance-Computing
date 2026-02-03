#!/usr/bin/env python3
"""
Script visualize kết quả LBM
Vẽ trường vận tốc và profile vận tốc
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

def load_lbm_data(filename):
    """Đọc file kết quả LBM"""
    data = np.loadtxt(filename, skiprows=1)
    x = data[:, 0]
    y = data[:, 1]
    rho = data[:, 2]
    ux = data[:, 3]
    uy = data[:, 4]
    
    # Tìm kích thước lưới
    nx = len(np.unique(x))
    ny = len(np.unique(y))
    
    # Reshape thành ma trận 2D
    X = x.reshape(nx, ny)
    Y = y.reshape(nx, ny)
    RHO = rho.reshape(nx, ny)
    UX = ux.reshape(nx, ny)
    UY = uy.reshape(nx, ny)
    
    return X, Y, RHO, UX, UY, nx, ny

def plot_velocity_field(X, Y, UX, UY, nx, ny):
    """Vẽ trường vận tốc"""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
    
    # Tính độ lớn vận tốc
    U_magnitude = np.sqrt(UX**2 + UY**2)
    
    # ===== Plot 1: Contour của độ lớn vận tốc =====
    levels = np.linspace(0, U_magnitude.max(), 20)
    cf = ax1.contourf(X, Y, U_magnitude, levels=levels, cmap='jet')
    ax1.set_xlabel('X', fontsize=12)
    ax1.set_ylabel('Y', fontsize=12)
    ax1.set_title('Độ lớn vận tốc |u|', fontsize=14, fontweight='bold')
    ax1.set_aspect('equal')
    cbar1 = plt.colorbar(cf, ax=ax1)
    cbar1.set_label('|u|', fontsize=11)
    
    # ===== Plot 2: Vector vận tốc (quiver) =====
    # Lấy mẫu thưa để vẽ cho dễ nhìn
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
    print("✅ Đã lưu trường vận tốc vào: velocity_field.png")

def plot_velocity_profile(UX, ny):
    """Vẽ profile vận tốc theo chiều Y"""
    # Lấy profile ở giữa lưới
    mid_x = UX.shape[0] // 2
    ux_profile = UX[mid_x, :]
    y_coords = np.arange(ny)
    
    # Vẽ profile
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    # ===== Subplot 1: Profile vận tốc =====
    ax1.plot(ux_profile, y_coords, 'b-o', linewidth=2, markersize=3, label='LBM (Mô phỏng)')
    ax1.set_xlabel('Vận tốc u_x', fontsize=12)
    ax1.set_ylabel('Y', fontsize=12)
    ax1.set_title('Profile vận tốc theo chiều Y', fontsize=14, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    
    # Fit parabola: u(y) = a*y^2 + b*y + c
    # Bỏ qua 2 điểm biên (y=0 và y=ny-1) để fit tốt hơn
    y_inner = y_coords[1:-1]
    ux_inner = ux_profile[1:-1]
    
    # Fit dạng parabola: u = u_max * (1 - ((y-y_center)/half_width)^2)
    y_center = (ny - 1) / 2.0
    coeffs = np.polyfit(y_inner - y_center, ux_inner, 2)  # a*y^2 + b*y + c
    u_fit = np.polyval(coeffs, y_coords - y_center)
    
    ax1.plot(u_fit, y_coords, 'r--', linewidth=2, label='Fit Parabolic')
    ax1.legend(fontsize=10)
    
    # Thêm thông tin
    u_max = ux_profile.max()
    u_max_theory = coeffs[2]  # hệ số c (đỉnh parabola)
    ax1.text(0.02, ny*0.85, f'u_max (LBM) = {u_max:.4f}\nu_max (fit) = {u_max_theory:.4f}', 
             fontsize=10, bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    # ===== Subplot 2: Sai số so với parabola =====
    error = ux_profile - u_fit
    ax2.plot(error * 1000, y_coords, 'g-o', linewidth=2, markersize=3)  # Nhân 1000 để dễ thấy
    ax2.axvline(x=0, color='k', linestyle='--', alpha=0.5)
    ax2.set_xlabel('Sai số (×10⁻³)', fontsize=12)
    ax2.set_ylabel('Y', fontsize=12)
    ax2.set_title('Sai số: LBM - Fit Parabolic', fontsize=14, fontweight='bold')
    ax2.grid(True, alpha=0.3)
    
    # Thêm thông tin sai số
    rmse = np.sqrt(np.mean(error[1:-1]**2))
    max_error = np.max(np.abs(error[1:-1]))
    ax2.text(0.02, ny*0.85, f'RMSE = {rmse:.6f}\nMax Error = {max_error:.6f}', 
             fontsize=10, bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.5),
             transform=ax2.get_yaxis_transform())
    
    plt.tight_layout()
    plt.savefig('velocity_profile.png', dpi=300, bbox_inches='tight')
    print("✅ Đã lưu profile vận tốc vào: velocity_profile.png")

def plot_density_field(X, Y, RHO):
    """Vẽ trường mật độ"""
    plt.figure(figsize=(10, 5))
    levels = np.linspace(RHO.min(), RHO.max(), 20)
    cf = plt.contourf(X, Y, RHO, levels=levels, cmap='viridis')
    plt.xlabel('X', fontsize=12)
    plt.ylabel('Y', fontsize=12)
    plt.title('Trường mật độ ρ', fontsize=14, fontweight='bold')
    plt.colorbar(cf, label='ρ')
    plt.axis('equal')
    plt.tight_layout()
    plt.savefig('density_field.png', dpi=300, bbox_inches='tight')
    print("✅ Đã lưu trường mật độ vào: density_field.png")

def main():
    """Hàm chính"""
    filename = 'lbm_result.dat'
    
    print(f"Đang đọc file: {filename}...")
    try:
        X, Y, RHO, UX, UY, nx, ny = load_lbm_data(filename)
        print(f"✅ Đã đọc lưới {nx}×{ny}")
        
        print("\nĐang vẽ biểu đồ...")
        plot_velocity_field(X, Y, UX, UY, nx, ny)
        plot_velocity_profile(UX, ny)
        plot_density_field(X, Y, RHO)
        
        # In thống kê
        print("\n" + "="*50)
        print("THỐNG KÊ KẾT QUẢ")
        print("="*50)
        print(f"Kích thước lưới: {nx} × {ny}")
        print(f"Mật độ trung bình: {RHO.mean():.6f}")
        print(f"Vận tốc u_x: min={UX.min():.6f}, max={UX.max():.6f}, avg={UX.mean():.6f}")
        print(f"Vận tốc u_y: min={UY.min():.6f}, max={UY.max():.6f}, avg={UY.mean():.6f}")
        print(f"Độ lớn vận tốc max: {np.sqrt(UX**2 + UY**2).max():.6f}")
        print("="*50)
        
    except FileNotFoundError:
        print(f"❌ Lỗi: Không tìm thấy file {filename}")
        print("   Hãy chạy ./lbm_serial trước để tạo file kết quả.")
    except Exception as e:
        print(f"❌ Lỗi: {e}")

if __name__ == '__main__':
    main()
