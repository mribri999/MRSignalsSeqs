## Rad229_MRI_Phantom – Create an object with MRI related properties.
#
# SYNTAX  - P, M = rad229_mri_phantom(acq)
#
# INPUTS  - acq is a structure that needs to contain:
#               acq['Nx'] pixels (number of columns)
#               acq['Ny'] pixels (number of rows)
#
# OUTPUTS - P - Is the phantom object matrix  that is acq['Nx'] by acq['Ny']
#           M - A 3D logical mask matrix. Each layer is an object in the phantom.
#
# EXAMPLE - acq = {'Nx': 128}
#           acq = {'Ny': 129}
#           P, M = rad229_mri_phantom(acq)
#
# DBE@STANFORD.EDU (April 2025) for Rad229

import numpy as np
from skimage.data import shepp_logan_phantom
from skimage.transform import resize

def rad229_mri_phantom(acq=None):
    if acq is None:
        acq = {'Nx': 128} # Define the matrix size (will be square)
        acq = {'Ny': 129}
    
    Nx = acq['Nx']
    Ny = acq['Ny']

    # Generate base phantom
    P = resize(shepp_logan_phantom(), (Ny, Nx), mode='reflect', anti_aliasing=True)

    # Phantom ellipses parameters (matching MATLAB's 'modified shepp-logan')
    # Parameters: [intensity, a, b, x0, y0, phi] similar to MATLAB phantom
    # Here we mimic the Shepp-Logan phantom components
    ellipses = [
        [1.0,   .69, .92,   0,    0,   0],
        [-.8,  .6624, .874,  0,  -.0184,   0],
        [-.2,  .1100, .310,  .22, 0, -18],
        [-.2,  .1600, .410, -.22, 0,  18],
        [0.1,  .2100, .250,  0,   .35,  0],
        [0.1,  .0460, .046,  0,   .1,   0],
        [0.1,  .0460, .046,  0,  -.1,   0],
        [0.1,  .0460, .023, -.08, -.605, 0],
        [0.1,  .0230, .023,  0,  -.606, 0],
        [0.1,  .0230, .046,  .06, -.605, 0]
    ]

    def ellipse_mask(params, shape):
        """Create binary mask for a single ellipse."""
        intensity, a, b, x0, y0, phi = params
        rows, cols = shape
        y, x = np.ogrid[-1:1:complex(rows), -1:1:complex(cols)]

        x_rot = (x - x0) * np.cos(np.deg2rad(phi)) + (y - y0) * np.sin(np.deg2rad(phi))
        y_rot = -(x - x0) * np.sin(np.deg2rad(phi)) + (y - y0) * np.cos(np.deg2rad(phi))

        return (x_rot / a)**2 + (y_rot / b)**2 <= 1

    # Create masks for each object
    #M = np.zeros((Ny, Nx, len(ellipses)), dtype=bool)
    #for i, e in enumerate(ellipses):
    #    mask = ellipse_mask([1] + e[1:], (Ny, Nx))  # Force intensity to 1
    #    M[:, :, i] = mask

    # Create masks for each object
    M = np.zeros((Ny, Nx, len(ellipses)), dtype=bool)
    for i, e in enumerate(ellipses):
        mask = ellipse_mask([1] + e[1:], (Ny, Nx))  # Force intensity to 1
        M[:, :, i] = np.flipud(mask)  # Flip up-down to match P

    # Parse the background and tissues
    Q = np.zeros((Ny, Nx, 3), dtype=bool)
    Q[:, :, 0] = ~M[:, :, 0] & ~M[:, :, 1]     # Background
    Q[:, :, 1] = M[:, :, 0] & M[:, :, 1]       # White-matter
    Q[:, :, 2] = M[:, :, 0] & ~M[:, :, 1]      # Skull

    # Final mask stack: [background, white-matter, skull, rest...]
    M_final = np.concatenate((Q, M[:, :, 2:]), axis=2)

    # Avoid overlapping in white matter
    overlap = np.sum(M_final, axis=2) > 1
    M_final[:, :, 1] = ~overlap & Q[:, :, 1]

    return P, M_final