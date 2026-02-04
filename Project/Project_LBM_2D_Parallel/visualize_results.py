#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

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
    
    levels = np.linspace(0, U_magnitude.max(), 20)
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
    print("✅ Đã lưu trường vận tốc vào: velocity_field.png")

def plot_velocity_profile(UX, ny):
    mid_x = UX.shape[0] // 2
    ux_profile = UX[mid_x, :]
    y_coords = np.arange(ny)
    
    plt.figure(figsize=(8, 6))
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
    plt.text(0.02, ny*0.85, f'u_max (LBM) = {u_max:.4f}', 
             fontsize=10, bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5),
             transform=plt.gca().get_yaxis_transform())
    
    plt.tight_layout()
    plt.savefig('velocity_profile.png', dpi=300, bbox_inches='tight')
    print("✅ Đã lưu profile vận tốc vào: velocity_profile.png")

def plot_density_field(X, Y, RHO):
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
    filename = 'lbm_result.dat'
    
    print(f"Đang đọc file: {filename}...")
    try:
        X, Y, RHO, UX, UY, nx, ny = load_lbm_data(filename)
        print(f"✅ Đã đọc lưới {nx}×{ny}")
        
        print("\nĐang vẽ biểu đồ...")
        plot_velocity_field(X, Y, UX, UY, nx, ny)
        plot_velocity_profile(UX, ny)
        plot_density_field(X, Y, RHO)
        
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