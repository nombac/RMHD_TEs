#!/usr/bin/env python3
"""
Python implementation of IDL's READU function for reading binary data files.
"""

import numpy as np
import os
import sys


def readu(filename, i=None, j=None, k=None, endian='big'):
    """
    Read binary data from a file with specified dimensions.
    
    Parameters:
    -----------
    filename : str
        Path to the binary file to read
    i : int, optional
        First dimension size
    j : int, optional
        Second dimension size  
    k : int, optional
        Third dimension size
    endian : str, optional
        Byte order ('big' or 'little'), default is 'big'
        
    Returns:
    --------
    numpy.ndarray
        Array containing the read data with specified shape
        Returns None if file cannot be opened
    """
    
    # Determine array shape based on provided dimensions
    if i is not None and j is not None and k is not None:
        shape = (i, j, k)
    elif i is not None and j is not None:
        shape = (i, j)
    elif i is not None:
        shape = (i,)
    else:
        print('% READU: <ERROR> KEYWORD IS WRONG')
        return None
    
    # Check if file exists
    if not os.path.exists(filename):
        print(f'% READU: <ERROR> COULD NOT OPEN {filename} (IGNORED)')
        return np.zeros(shape, dtype=np.float32)
    
    try:
        # Determine byte order
        byte_order = '>' if endian == 'big' else '<'
        
        # Read binary data
        with open(filename, 'rb') as f:
            # Skip first 4 bytes (data_start=4 in IDL)
            f.seek(4)
            
            # Read float32 data
            dtype = np.dtype(f'{byte_order}f4')
            data = np.fromfile(f, dtype=dtype, count=np.prod(shape))
            
            # Reshape to specified dimensions
            # Note: IDL uses column-major (Fortran) order, NumPy uses row-major (C) order
            # So we need to reshape with Fortran order and potentially transpose
            array = data.reshape(shape, order='F')
            
        return array
        
    except Exception as e:
        print(f'% READU: <ERROR> COULD NOT READ {filename}: {e}')
        return np.zeros(shape, dtype=np.float32)


def main():
    """
    Test the readu function with vave.data file.
    """
    
    # Check if test/h/vave.data exists
    vave_path = 'test/h/vave.data'
    resolv_path = 'test/h/resolv.data'
    
    if os.path.exists(vave_path):
        print(f"Found vave.data at {vave_path}")
        
        # First, try to read resolv.data to get dimensions
        n_t = None
        n_vave = None
        
        if os.path.exists(resolv_path):
            print(f"Reading dimensions from {resolv_path}")
            try:
                with open(resolv_path, 'r') as f:
                    line = f.readline().strip()
                    parts = line.split()
                    if len(parts) >= 2:
                        n_t = int(parts[0])
                        n_vave = int(parts[1])
                        print(f"Dimensions: n_t={n_t}, n_vave={n_vave}")
            except Exception as e:
                print(f"Could not read resolv.data: {e}")
        
        # If we have dimensions, try to read vave.data
        if n_t and n_vave:
            print(f"\n--- Reading vave.data as {n_t}x{n_vave} array ---")
            data = readu(vave_path, i=n_t, j=n_vave)
            if data is not None:
                print(f"Successfully read vave.data")
                print(f"Shape: {data.shape}")
                print(f"Data type: {data.dtype}")
                print(f"Min value: {np.min(data)}")
                print(f"Max value: {np.max(data)}")
                print(f"Mean value: {np.mean(data)}")
                print(f"\nFirst 5 time steps (first 5 rows):")
                print(data[:5, :min(10, n_vave)])
                
                # Check for NaN or Inf values
                nan_count = np.sum(np.isnan(data))
                inf_count = np.sum(np.isinf(data))
                print(f"\nNaN values: {nan_count}")
                print(f"Inf values: {inf_count}")
                
                # Save a sample for verification
                print("\n--- Saving sample data for verification ---")
                np.save('vave_data_sample.npy', data)
                print("Saved data to vave_data_sample.npy")
        else:
            print("Warning: Could not determine dimensions from resolv.data")
            print("Attempting to read with guessed dimensions...")
            
            # Try to get file size to estimate dimensions
            file_size = os.path.getsize(vave_path)
            # Subtract 4 bytes for header, divide by 4 bytes per float
            n_floats = (file_size - 4) // 4
            print(f"File size: {file_size} bytes")
            print(f"Estimated number of float values: {n_floats}")
            
            # Try reading as 1D array first
            print(f"\n--- Reading as 1D array ({n_floats} elements) ---")
            data_1d = readu(vave_path, i=n_floats)
            if data_1d is not None:
                print(f"Shape: {data_1d.shape}")
                print(f"First 20 values: {data_1d[:20]}")
                print(f"Min: {np.min(data_1d)}, Max: {np.max(data_1d)}")
    else:
        print(f"vave.data not found at {vave_path}")
        print("Please place vave.data file in test/h/ directory")
        print("\nRunning basic functionality tests instead...")
        
        # Create test binary file
        test_filename = 'test_data.bin'
        test_data = np.arange(24, dtype=np.float32)
        
        with open(test_filename, 'wb') as f:
            f.write(b'\x00\x00\x00\x00')  # 4-byte header
            test_data.astype('>f4').tofile(f)  # Big-endian float32
        
        print("\n--- Test: Reading test data ---")
        result = readu(test_filename, i=24)
        if result is not None:
            print(f"Shape: {result.shape}")
            print(f"Data matches: {np.allclose(result, test_data)}")
        
        # Cleanup
        if os.path.exists(test_filename):
            os.remove(test_filename)
    
    print("\n✓ Test completed")


if __name__ == "__main__":
    main()