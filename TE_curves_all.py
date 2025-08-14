#!/usr/bin/env python3
"""
S-curve plotting script for Zeus radiation hydrodynamics simulations.
Converted from IDL s_curve_all.pro to Python.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import os
from target_config import TARGETS, COLORS, DATA_BASE_PATH
from zeus_param import zeus_param
from readu import readu

# Physical constants (matching IDL startup.pro)
SIGMAB = 5.67040e-5  # Stefan-Boltzmann constant
YEAR = 3.1536e7      # Year in seconds
MSOL = 1.989e33      # Solar mass


def main():
    """Main function to process targets and create plots."""
    
    print("Processing S-curve targets...")
    
    # Storage for plot data
    plot_data = []
    
    # Process each target
    for j, targ in enumerate(TARGETS):
        # Skip dummy entry
        if targ['name'] == 'dummy':
            break
            
        print(f"Processing target: {targ['name']}")
        
        # Initialize variables
        sigma, teff, omega, lz, pres, stress = 0, 0, 0, 0, 0, 0
        
        # Handle special cases (h2006, h2007) with hardcoded values
        if targ['name'] == 'h2006': # DOI 10.1086/499153
            omega = 5.9  # 回転角速度 ω=5.90 s^-1（p.3 Table 1）
            sigma = 9.89e4  # 表面密度 Σ=9.89×10^4 g cm^-2（p.3 Table 1）
            teff = 5.3e5  # 有効温度 5.3×10^5 K（p.5：time-averaged effective temperature）
            flux = SIGMAB * (teff)**4  # ステファン＝ボルツマン則から放射フラックス計算
            alpha = 0.016  # 時間平均αパラメータ（p.6）
            stress = 2.0/3.0 * flux / omega  # 放射フラックスとωから算出したストレス
            pres = stress / alpha  # ストレスとαから全圧を計算
            lz = 1.0  # 高さ方向スケーリング因子（仮定値）
        elif targ['name'] == 'h2007': # DOI 10.1086/519515
            omega = 17.0  # 回転角速度 ω=17 s^-1（p.2）
            sigma = 4.7e4  # 表面密度 Σ=4.7×10^4 g cm^-2（p.2）
            flux0 = SIGMAB * (9e5)**4  # 有効温度9.0×10^5 Kからの放射フラックス（p.2）
            flux = 1.65 * flux0  # 初期値より65%大きい放射フラックス（p.4）
            teff = (flux / 2.0 / SIGMAB)**0.25  # 放射フラックスから片面の有効温度を計算
            alpha = 0.03  # 時間平均αパラメータ（p.4）
            stress = 2.0/3.0 * flux / omega  # 放射フラックスとωから算出したストレス
            pres = stress / alpha  # ストレスとαから全圧を計算
            lz = 1.0  # 高さ方向スケーリング因子（仮定値）
            
        else:
            # Regular targets - read from simulation data
            try:
                # Construct path
                path = os.path.join(DATA_BASE_PATH, targ["name"])
                
                # Read simulation parameters
                omega, h0, gamma, mmw, isoth, torb, lx, ly, lz = zeus_param(path)
                
                # Read volume-averaged history data
                resolv_file = os.path.join(path, 'h', 'resolv.data')
                with open(resolv_file, 'r') as f:
                    line = f.readline().strip().split()
                    n_t = int(line[0])
                    n_vave = int(line[1])
                
                # Read vave.data
                vave_file = os.path.join(path, 'h', 'vave.data')
                a = readu(vave_file, i=n_t, j=n_vave)
                
                # Calculate time averaging window
                if targ['ave'][0] == 0 and targ['ave'][1] != 0:
                    nmin = n_t - 1 - targ['ave'][1] * 100
                    nmax = n_t - 1
                else:
                    nmin = int(targ['ave'][0] / 0.01)
                    nmax = int(targ['ave'][1] / 0.01)
                
                # Apply bounds
                if nmin < 0:
                    nmin = 0
                if nmax > n_t - 1:
                    nmax = n_t - 1
                
                # Calculate thermal pressure
                gammar = 4.0 / 3.0
                prad = (gammar - 1.0) * a[:, 10] * lz  # 11-1 in IDL (1-based) = 10 in Python (0-based)
                
                # Calculate gas pressure based on target type
                color_key = targ['color_key']
                if color_key == 'h2009':
                    gamma = 5.0 / 3.0
                    pgas = (gamma - 1.0) * a[:, 3] * lz  # 4-1 in IDL = 3 in Python
                elif color_key == 'h2011':
                    gamma = 1.4
                    pgas = (gamma - 1.0) * a[:, 3] * lz
                elif color_key in ['h2014', 'h2015', 'h2016']:
                    pgas = a[:, 80] * lz  # 81-1 in IDL = 80 in Python
                else:
                    # Default case
                    pgas = (gamma - 1.0) * a[:, 3] * lz
                
                pres = prad + pgas
                pres = np.mean(pres[nmin:nmax])
                
                # Calculate surface density
                sigma = a[:, 2] * lz  # 3-1 in IDL = 2 in Python
                sigma = np.mean(sigma[nmin:nmax])
                
                # Calculate flux from shear stress work
                if color_key == 'h2016':
                    # Self-gravity case
                    flux = (a[:, 88] +    # Maxwell stress (89-1 in IDL = 88 in Python)
                           a[:, 87] +     # Reynolds stress (88-1 in IDL = 87 in Python) 
                           a[:, 78])      # Self-gravity (79-1 in IDL = 78 in Python)
                    flux *= lz
                else:
                    # Standard case
                    flux = 1.5 * omega * (a[:, 11] + a[:, 12]) * lz  # 12-1, 13-1 in IDL = 11, 12 in Python
                
                flux = np.mean(flux[nmin:nmax])
                teff = (flux / 2.0 / SIGMAB)**0.25
                
                # Calculate stress for alpha plot
                stress = 2.0/3.0 * flux / omega
                
            except Exception as e:
                print(f"Error processing {targ['name']}: {e}")
                continue
        
        # Store data for plotting (only if values are valid)
        if sigma > 0 and teff > 0:
            plot_data.append({
                'name': targ['name'],
                'sigma': sigma,
                'teff': teff,
                'pres': pres,
                'stress': stress,
                'omega': omega,
                'lz': lz,
                'color_key': targ['color_key']
            })
        
        # Output the results (matching IDL print statement)
        print(f"{targ['name']}: sigma={sigma:.2e}, teff={teff:.2e}, omega={omega:.2e}, lz={lz:.2e}")
    
    # Create plots
    create_plots(plot_data)


def create_plots(plot_data):
    """Create both S-curve and alpha plots."""
    
    # Create S-curve plot
    print("\nCreating S-curve plot...")
    fig1, ax1 = create_scurve_plot(plot_data)
    
    # Create alpha plot  
    print("Creating alpha plot...")
    fig2, ax2 = create_alpha_plot(plot_data)


def create_scurve_plot(plot_data):
    """Create S-curve plot (sigma vs teff)."""
    
    # Create output directory if it doesn't exist
    os.makedirs('./outputs', exist_ok=True)
    
    # Plot settings
    output_file = './outputs/TE_curves_all.pdf'
    xlabel = r'Surface density [g cm$^{-2}$]'
    ylabel = r'Effective temperature [K]'
    xrange = [1e1, 1e6]
    yrange = [1e1, 1e7]
    
    # Create figure with IDL-like proportions
    xsize = 8.5  # inches
    dxw = 0.75
    dyw = 0.85
    ysize = xsize * dxw / np.sqrt(2.0) / dyw * 1.5
    
    fig, ax = plt.subplots(figsize=(xsize, ysize))
    
    # Set up axes
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(xrange)
    ax.set_ylim(yrange)
    ax.set_xlabel(xlabel, fontsize=14)
    ax.set_ylabel(ylabel, fontsize=14)
    ax.tick_params(axis='both', which='major', labelsize=12)
    ax.tick_params(axis='both', which='minor', labelsize=10)
    
    # Plot data points
    for data in plot_data:
        color = COLORS.get(data['color_key'], 'black')
        ax.scatter(data['sigma'], data['teff'], 
                  c=color, s=50, marker='o', alpha=0.8, 
                  edgecolors='black', linewidths=0.5,
                  label=data['color_key'] if data == plot_data[0] else "")
    
    # Save and show
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"S-curve plot saved to: {output_file}")
    plt.show()
    
    return fig, ax


def create_alpha_plot(plot_data):
    """Create alpha plot (pres vs stress)."""
    
    # Create output directory if it doesn't exist
    os.makedirs('./outputs', exist_ok=True)
    
    # Plot settings
    output_file = './outputs/alpha_all.pdf'
    xlabel = r'Pressure [dyn cm$^{-2}$]'
    ylabel = r'Stress [dyn cm$^{-2}$]'
    
    # Determine ranges from data
    valid_pres = [d['pres'] for d in plot_data if d['pres'] > 0]
    valid_stress = [d['stress'] for d in plot_data if d['stress'] > 0]
    
    if valid_pres and valid_stress:
        xrange = [min(valid_pres) * 0.1, max(valid_pres) * 10]
        yrange = [min(valid_stress) * 0.1, max(valid_stress) * 10]
    else:
        xrange = [1e9, 1e23]
        yrange = [1e9, 1e23]
    
    # Create figure
    xsize = 8.5
    dxw = 0.75
    dyw = 0.85
    ysize = xsize * dxw / np.sqrt(2.0) / dyw * 1.5
    
    fig, ax = plt.subplots(figsize=(xsize, ysize))
    
    # Set up axes
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(xrange)
    ax.set_ylim(yrange)
    ax.set_xlabel(xlabel, fontsize=14)
    ax.set_ylabel(ylabel, fontsize=14)
    ax.tick_params(axis='both', which='major', labelsize=12)
    ax.tick_params(axis='both', which='minor', labelsize=10)
    
    # Plot reference lines (alpha = const) with different line styles
    xvals = np.logspace(np.log10(xrange[0]), np.log10(xrange[1]), 100)
    ax.plot(xvals, xvals * 1e-1, 'k-', alpha=0.5, linewidth=1, label=r'$\alpha=0.1$')      # solid line
    ax.plot(xvals, xvals * 1e-2, 'k--', alpha=0.5, linewidth=1, label=r'$\alpha=0.01$')    # dashed line
    ax.plot(xvals, xvals * 1e-3, 'k:', alpha=0.5, linewidth=1, label=r'$\alpha=0.001$')    # dotted line
    
    # Plot data points
    for data in plot_data:
        if data['pres'] > 0 and data['stress'] > 0:
            color = COLORS.get(data['color_key'], 'gray')
            # Plot gray background point first (larger)
            ax.scatter(data['pres'], data['stress'], c='gray', s=60, marker='o', alpha=0.6)
            # Plot colored point on top (smaller)
            ax.scatter(data['pres'], data['stress'], c=color, s=40, marker='o', alpha=0.8)
    
    # Only add legend if there are reference lines
    if len(ax.lines) > 0:
        ax.legend(loc='upper left', fontsize=10)
    
    # Save and show
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Alpha plot saved to: {output_file}")
    plt.show()
    
    return fig, ax


if __name__ == '__main__':
    main()
