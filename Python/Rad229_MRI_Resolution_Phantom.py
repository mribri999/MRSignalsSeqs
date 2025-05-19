# This function creates a T2* resolution phantom with a range of feature sizes.
# The input matrix size is flexible, but the number of objects is fixed at five. 
# The resolution elements have single, double, and quadruple line spacing.
# 
# DBE@STANFORD.EDU (April 2025) Python version for Rad229

import numpy as np

def Rad229_MRI_Resolution_Phantom(acq=None, obj=None):
    # Default values if no input is provided
    if acq is None or obj is None:
        acq = {}  # Empty dictionary for acq
        acq['Nx'] = 128  # Define Nx (number of columns)
        acq['Ny'] = 256  # Define Ny (number of rows)
        obj = {'T2star': [10e-3, 25e-3, 50e-3, 100e-3, 1000e-3]}

    n_obj = len(obj['T2star'])  # Number of objects [hard coded]
   # n_elem = 9 # Number of elements (lines) per object

    gap = int(acq['Nx'] / 16)  # Gap between objects and edges
    wid = int(3 * gap)  # Width (number of columns) for objects
    P = np.zeros((acq['Ny'], acq['Nx']))  # Initialize an array

    # Define number of elements to fill the y-direction
    n_elem = int((acq['Ny'] - 4 * gap)/14) + 1  # 4 gaps for 3 resolution elements that 2+4+8=14 pixels each
    
    # Subroutine to define the phantom values
    def define_obj_values(P_obj, obj, n_obj, n_elem, r0, r1, c0, c1):
        for c in range(n_obj):
            for r in range(n_elem):
                P_obj[r0[r]:r1[r], c0[c]:c1[c]] = obj['T2star'][c]
        return P_obj

    # Define the columne spacing (fixed for all objects)
    c0 = np.arange(gap + 1, acq['Nx'], wid)  # Start column index
    c1 = c0 + 2 * gap   # End column index # NOT SURE WHY THIS IS 2*gap...

    # Define the fine line (row) spacing
    r0 = gap + np.arange(0, 2 * n_elem, 2)  # Start row index
    r1 = r0 + 1   # End row index
    P = define_obj_values(P, obj, n_obj, n_elem, r0, r1, c0, c1)

    # Define the medium line (row) spacing
    r0 = np.max(r1) + gap + np.arange(0, 4 * n_elem, 4)  # Start row index
    r1 = r0 + 2   # End row index
    P = define_obj_values(P, obj, n_obj, n_elem, r0, r1, c0, c1)

    # Define the coarse line (row) spacing
    r0 = np.max(r1) + gap + np.arange(0, 8 * n_elem, 8)  # Start row index
    r1 = r0 + 4   # End row index
    P = define_obj_values(P, obj, n_obj, n_elem, r0, r1, c0, c1)

    return P