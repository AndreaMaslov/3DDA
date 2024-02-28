"""
Droplet density analysis from LAMMPS dump data - Angels between ellipsoid tangent plane and reference horizontal plane in intersection points module
"""

# Standart python packages
import math
import numpy as np

def ellipsoid_tangent_angle(frame_number, el_center, el_radii, x_intersection_point, y_intersection_point, z_intersection_point):

    # Proved fact:
    #   The angle between two 3D planes is same as the angle between their normal vectors
    #
    # Consider:
    #   Three unit vectors: i(1,0,0), j(0,1,0) and k(0,0,1)
    #   Fitted ellipsoid surface: F(x,y,z) = 0
    #   Point P(x0,y0,z0) is intersection point of z_ref plane and fitted ellipsoid surface
    #
    # Then:
    #   Vector k is normal to the z_ref plane at any point, ergo also at point P(x0,y0,z0)
    #   Gradient vector of F(x,y,z) = dF/dx*i + dF/dy*j + dF/dz*k at point P is normal to the ellipsoid tangent plane at point P 
    #   Angle between two normal vectors "k" and "dF/dx*i + dF/dy*j + dF/dz*k" at point P (||n|| --> magnitude of vector n):
    #      cos(angle)=dot_product(n1,n2)/(||n1||*||n2||) where n1 is "dF/dx*i + dF/dy*j + dF/dz*k" and n2 is "k" with ||"k"||=1
    #      angle=arccos(dot_product("dF/dx*i + dF/dy*j + dF/dz*k","k")/||"dF/dx*i + dF/dy*j + dF/dz*k"||)

    # Ellipsoid equation F(x,y,z) = 0:
    #   (x - el_center[0])^2/el_radii[0]^2 + (y - el_center[1])^2/el_radii[1]^2 + (z - el_center[2])^2/el_radii[2]^2 - 1 = 0
    # Gradient vector:
    #   (2*(x - el_center[0])/el_radii[0]^2)*i + (2*(y - el_center[1])/el_radii[1]^2)*j + (2*(z - el_center[2])/el_radii[2]^2)*k

    # Calculate angle between two normal vectors for every intersection point
    lst_tangent_angles_per_frame = []
    for i in range(z_intersection_point.shape[0]):
        #
        magnitude_gradient_vector = math.sqrt((2*(x_intersection_point[i] - el_center[0])/el_radii[0]**2)**2 +
                                              (2*(y_intersection_point[i] - el_center[1])/el_radii[1]**2)**2 +
                                              (2*(z_intersection_point[i] - el_center[2])/el_radii[2]**2)**2)
        #
        dot_product_gradient_vector_k_vector = 2*(z_intersection_point[i] - el_center[2])/el_radii[2]**2
        #print("dot_product_gradient_vector_k_vector: ", dot_product_gradient_vector_k_vector)
        #
        angle_in_radians = math.acos(dot_product_gradient_vector_k_vector/magnitude_gradient_vector)  #  domain of acos result 0 <------> pi
        #print("Tangent angle in radians: ", angle_in_radians)
        #
        angle_in_degrees = math.degrees(angle_in_radians)
        #print("Tangent angle in degrees: ", angle_in_degrees)

        # Add angle to the list
        lst_tangent_angles_per_frame.append(angle_in_degrees)

    #print("Frame: ", frame_number)
    #print("Number of angles: ", len(lst_angles_per_frame))
    #print("Angles values: ", lst_tangent_angles_per_frame[:])

    return lst_tangent_angles_per_frame

