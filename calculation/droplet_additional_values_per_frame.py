"""
Droplet density analysis from LAMMPS dump data -Droplet additional values per frame module
"""

# Standart python packages
import numpy as np
import math


def droplet_additional_values(frame_number, z_ref_ellipse_params, el_center, el_radii, z_ref, lst_radius_values, lst_area_values, lst_volume_values):

    # Calculate Droplet Additional Ellipse Values
    radius_value, area_value = get_droplet_additional_ellipse_values(z_ref_ellipse_params)
    #print("radius_value: ", radius_value)
    #print("area_value: ", area_value)
    lst_radius_values.append(radius_value)
    lst_area_values.append(area_value)

    # Calculate Droplet Additional Ellipsoid Values
    volume_value = get_droplet_additional_ellipsoid_values(el_center, el_radii, z_ref)
    #print("volume_value: ", volume_value)
    lst_volume_values.append(volume_value)

    return lst_radius_values, lst_area_values, lst_volume_values

def get_droplet_additional_ellipse_values(z_ref_ellipse_params):

    x0, y0, ap, bp, e, phi = z_ref_ellipse_params
    # ap - ellipse semi-major axis
    # bp - ellipse semi-minor axis
    #print("ap: ", ap)
    #print("bp: ", bp)

    radius_value = (ap + bp) / 2.0
    area_value = math.pi * ap * bp
    
    return radius_value, area_value

def get_droplet_additional_ellipsoid_values(el_center, el_radii, z_ref):
    """
    Based on: https://keisan.casio.com/exec/system/1311572253

    """
    a = el_radii[0]
    b = el_radii[1]
    c = el_radii[2]
    z0 = el_center[2]
    #
    #print("a: ", a)
    #print("b: ", b)
    #print("c: ", c)
    #print("z0: ", z0)
    #print("z_ref: ", z_ref)

    z_max = z0 + c  # top of ellipsoid
    h = z_max - z_ref  # ≦2c
    #
    #print("z_max: ", z_max)
    #print("h: ", h)

    volume_value = ((math.pi * a * b)/(3 * c**2)) * (h**2 * (3 * c - h))
    
    return volume_value