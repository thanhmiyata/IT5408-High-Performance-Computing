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
    plt.figure(figsize=(8, 6))
    plt.plot(ux_profile, y_coords, 'b-o', linewidth=2, markersize=4)
    plt.xlabel('Vận tốc u_x', fontsize=12)
    plt.ylabel('Y', fontsize=12)
    plt.title('Profile vận tốc theo chiều Y (tại giữa lưới)', 
             fontsize=14, fontweight='bold')
    plt.grid(True, alpha=0.3)
    
    # Thêm đường tham chiếu (profile Poiseuille lý thuyết)
    # u(y) = u_max * (1 - (2y/H - 1)^2) cho dòng chảy Poiseuille
    y_norm = (y_coords - ny/2) / (ny/2)
    u_poiseuille = ux_profile.max() * (1 - y_norm**2)
    plt.plot(u_poiseuille, y_coords, 'r--', linewidth=2, label='Lý thuyết (Poiseuille)')
    plt.legend(fontsize=10)
    
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
