"""
Droplet density analysis from LAMMPS dump data - Convex Hull of density points
"""

# Standart python packages
import numpy as np

# Packages to install with "pip3 install -U package-name"
#
#   Transformation
#
from scipy.spatial import ConvexHull


def density_convex_hull(dens_nonzero_xyz, z_ref, z_com):

    #print("Statistics for dens_nonzero_xyz ndarray:")
    #print("    Number of axes (dimensions): ", dens_nonzero_xyz.ndim)
    #print("    Number of elements in each dimension: ", dens_nonzero_xyz.shape)
    #print("    All dens_nonzero_xyz elements: ", dens_nonzero_xyz)

    #    Perform Convex hull of rz_coords points
    #
    hull = ConvexHull(dens_nonzero_xyz)
    #
    #    Get ndarray of points (xyz coordinates) in the convex hull.
    hull_vertices_indx = hull.vertices
    #print("Statistics for hull_vertices_indx ndarray:")
    #print("    Number of axes (dimensions): ", hull_vertices_indx.ndim)
    #print("    Number of elements in each dimension: ", hull_vertices_indx.shape)
    #print("    All hull_vertices_indx elements: ", hull_vertices_indx)
    hull_vertices_xyz_full = hull.points[hull.vertices]
    #print("Statistics for hull_vertices_xyz_full ndarray:")
    #print("    Number of axes (dimensions): ", hull_vertices_xyz_full.ndim)
    #print("    Number of elements in each dimension: ", hull_vertices_xyz_full.shape)
    #print("    All hull_vertices_xyz_full elements: ", hull_vertices_xyz_full)
    #
    #   Consider only hull vertices above z_ref plane and below z_com
    hull_vertices_xyz = hull_vertices_xyz_full[(hull_vertices_xyz_full[:, 2] >= z_ref) & (hull_vertices_xyz_full[:, 2] <= z_com)]
    #print("Statistics for hull_vertices_xyz ndarray:")
    #print("    Number of axes (dimensions): ", hull_vertices_xyz.ndim)
    #print("    Number of elements in each dimension: ", hull_vertices_xyz.shape)
    #print("    All hull_vertices_xyz elements: ", hull_vertices_xyz)
    #
    #    Get X, Y coordinats from ndarray of hull_vertices_xyz
    #hull_vertices_x = hull_vertices_xyz[:, 0]
    #hull_vertices_y = hull_vertices_xyz[:, 1]
    #hull_vertices_z = hull_vertices_xyz[:, 2]

    return hull_vertices_indx, hull_vertices_xyz
