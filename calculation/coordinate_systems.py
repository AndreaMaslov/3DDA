"""
Droplet density analysis from LAMMPS dump data - Coordinate systems module
"""

# Standart python packages
import numpy as np

def z_ref(dens_nonzero_xyz, args_zv):

    #   Get reference horizontal plane
    #
    z_min_dens_nonzero_xyz = np.min(dens_nonzero_xyz[:, 2])
    #print("Min value of all dens_nonzero_xyz z coordinates: ", z_min_dens_nonzero_xyz)
    ref_z = z_min_dens_nonzero_xyz + args_zv
    #print("Reference horizontal plane:", ref_z)

    return ref_z
