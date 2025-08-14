#!/usr/bin/env python3
"""
Python implementation of IDL's ZEUS_PARAM procedure for reading Zeus simulation parameters.
"""

import os
import subprocess
import numpy as np


def spawnx(command, dtype='float'):
    """
    Execute a shell command and return the result as specified type.
    
    Parameters:
    -----------
    command : str
        Shell command to execute
    dtype : str
        Data type to convert result to ('long', 'float', 'double')
        
    Returns:
    --------
    Converted result from command output
    """
    try:
        result = subprocess.run(command, shell=True, capture_output=True, text=True)
        output = result.stdout.strip()
        
        if dtype == 'long':
            return int(output) if output else 0
        elif dtype in ['float', 'double']:
            return float(output) if output else 0.0
        else:
            return output
    except:
        return 0 if dtype == 'long' else 0.0


def zeus_param(homepath):
    """
    Read Zeus simulation parameters from configuration files.
    
    Parameters:
    -----------
    homepath : str
        Path to the Zeus simulation directory
        
    Returns:
    --------
    tuple : (omega, h0, gamma, mmw, isoth, torb, lx, ly, lz)
        omega : float - Angular frequency
        h0 : float - Scale height
        gamma : float - Adiabatic index
        mmw : float - Mean molecular weight
        isoth : int - Isothermal flag
        torb : float - Orbital period
        lx, ly, lz : float - Box dimensions
    """
    
    # PATHS
    initpath = os.path.join(homepath, 'init')
    
    # CHECK ISOTHERMAL/RADIATION
    isoth = spawnx(f"cat {homepath}/src/Makefile | grep 'PROBLEM =' | grep -v '#' | grep iso | wc -l", 'long')
    iokuz = spawnx(f"cat {homepath}/src/Makefile | grep 'PROBLEM =' | grep -v '#' | grep okuzumi | wc -l", 'long')
    iadia = spawnx(f"cat {homepath}/src/Makefile | grep 'PROBLEM =' | grep -v '#' | grep ad | wc -l", 'long')
    isano = spawnx(f"cat {homepath}/z3dinput | grep sanocon | wc -l", 'long')
    isoth += iokuz
    
    # SHEARING BOX PARAMETERS (OMEGA, H0)
    if isoth:
        omega = 1.0
        h0 = 1.0
    else:
        # Try to read from various parameter files
        param_files = ['iparax.data', 'normalize.data', 'ipara.data']
        omega = 1.0
        h0 = 1.0
        
        for param_file in param_files:
            filepath = os.path.join(initpath, param_file)
            if os.path.exists(filepath):
                try:
                    with open(filepath, 'r') as f:
                        line = f.readline().strip()
                        parts = line.split()
                        if len(parts) >= 2:
                            omega = float(parts[0])
                            h0 = float(parts[1])
                            break
                except:
                    continue
    
    # GAS PARAMETERS (GAMMA, MMW)
    if isoth:
        gamma = 2.0
        mmw = 1.0
    else:
        gamma = spawnx(f"cat {homepath}/z3dinput | grep gascon | awk '{{print $8}}'", 'float')
        mmw = spawnx(f"cat {homepath}/z3dinput | grep gascon | awk '{{print $4}}'", 'float')
    
    if isano and iadia:
        gamma = spawnx(f"cat {homepath}/z3dinput | grep sanocon | awk '{{print $8}}'", 'float')
    
    if isoth or isano:
        h00 = 1.0
    else:
        h00 = h0
    
    # BOX SIZE
    lx_max = spawnx(f"cat {homepath}/z3dinput | grep ggen1 | awk '{{print $12}}'", 'float')
    lx_min = spawnx(f"cat {homepath}/z3dinput | grep ggen1 | awk '{{print $8}}'", 'float')
    lx = (lx_max - lx_min) * h00
    
    ly_max = spawnx(f"cat {homepath}/z3dinput | grep ggen2 | awk '{{print $12}}'", 'float')
    ly_min = spawnx(f"cat {homepath}/z3dinput | grep ggen2 | awk '{{print $8}}'", 'float')
    ly = (ly_max - ly_min) * h00
    
    lz_max = spawnx(f"cat {homepath}/z3dinput | grep ggen3 | awk '{{print $12}}'", 'float')
    lz_min = spawnx(f"cat {homepath}/z3dinput | grep ggen3 | awk '{{print $8}}'", 'float')
    lz = (lz_max - lz_min) * h00
    
    # ORBITAL PERIOD
    torb = (2.0 * np.pi) / omega
    
    return omega, h0, gamma, mmw, isoth, torb, lx, ly, lz


def main():
    """
    Test the zeus_param function with actual Zeus files in test/ directory.
    """
    
    # Import readu function from readu.py
    from readu import readu
    
    # Use the test directory with actual Zeus files
    test_dir = 'test'
    
    # Check if test directory exists
    if not os.path.exists(test_dir):
        print(f"Error: {test_dir} directory not found")
        print("Please create the following structure:")
        print("  test/")
        print("  test/src/Makefile")
        print("  test/init/")
        print("  test/z3dinput")
        return
    
    # Check for required files
    required_files = [
        f"{test_dir}/src/Makefile",
        f"{test_dir}/z3dinput"
    ]
    
    for filepath in required_files:
        if os.path.exists(filepath):
            print(f"✓ Found: {filepath}")
        else:
            print(f"✗ Missing: {filepath}")
    
    # Check for optional parameter files
    param_files = [
        f"{test_dir}/init/iparax.data",
        f"{test_dir}/init/normalize.data",
        f"{test_dir}/init/ipara.data"
    ]
    
    print("\nOptional parameter files:")
    for filepath in param_files:
        if os.path.exists(filepath):
            print(f"✓ Found: {filepath}")
        else:
            print(f"  Not found: {filepath}")
    
    print("\n--- Reading Zeus parameters from test/ directory ---")
    try:
        omega, h0, gamma, mmw, isoth, torb, lx, ly, lz = zeus_param(test_dir)
        
        print(f"\nParameters extracted:")
        print(f"  omega (angular frequency) = {omega}")
        print(f"  h0 (scale height) = {h0}")
        print(f"  gamma (adiabatic index) = {gamma}")
        print(f"  mmw (mean molecular weight) = {mmw}")
        print(f"  isoth (isothermal flag) = {isoth}")
        print(f"  torb (orbital period) = {torb:.6f}")
        print(f"  lx (box size x) = {lx}")
        print(f"  ly (box size y) = {ly}")
        print(f"  lz (box size z) = {lz}")
        
        # Additional diagnostics
        print(f"\nDerived values:")
        print(f"  Box volume = {lx * ly * lz:.6f}")
        print(f"  Aspect ratio (lz/lx) = {lz/lx:.3f}")
        
        if isoth:
            print("\n  Running in ISOTHERMAL mode")
        else:
            print("\n  Running in RADIATION mode")
        
        # Try to calculate sigma from vave.data if available
        vave_path = f"{test_dir}/h/vave.data"
        resolv_path = f"{test_dir}/h/resolv.data"
        
        if os.path.exists(vave_path) and os.path.exists(resolv_path):
            print("\n--- Calculating surface density (sigma) ---")
            
            # Read dimensions from resolv.data
            try:
                with open(resolv_path, 'r') as f:
                    line = f.readline().strip()
                    parts = line.split()
                    if len(parts) >= 2:
                        n_t = int(parts[0])
                        n_vave = int(parts[1])
                        print(f"Time steps: {n_t}, Variables: {n_vave}")
                        
                        # Read vave.data
                        data = readu(vave_path, i=n_t, j=n_vave)
                        
                        if data is not None:
                            # Calculate sigma = a[*,3-1] * lz (column 3 in IDL = index 2 in Python)
                            sigma_array = data[:, 2] * lz  # Column 3 (0-indexed: column 2)
                            
                            # Calculate time average (using last 100 time steps or all if less)
                            n_avg = min(100, n_t)
                            sigma_mean = np.mean(sigma_array[-n_avg:])
                            sigma_std = np.std(sigma_array[-n_avg:])
                            
                            print(f"\nSurface density (sigma):")
                            print(f"  Mean (last {n_avg} steps): {sigma_mean:.6e}")
                            print(f"  Std deviation: {sigma_std:.6e}")
                            print(f"  Min: {np.min(sigma_array):.6e}")
                            print(f"  Max: {np.max(sigma_array):.6e}")
                            
                            # Show evolution
                            print(f"\n  First 5 time steps: {sigma_array[:5]}")
                            print(f"  Last 5 time steps: {sigma_array[-5:]}")
                            
            except Exception as e:
                print(f"Error calculating sigma: {e}")
        else:
            print(f"\nNote: vave.data or resolv.data not found in {test_dir}/h/")
            print("      Cannot calculate surface density (sigma)")
        
    except Exception as e:
        print(f"Error reading parameters: {e}")
        import traceback
        traceback.print_exc()
    
    print("\n✓ Test completed")


if __name__ == "__main__":
    main()