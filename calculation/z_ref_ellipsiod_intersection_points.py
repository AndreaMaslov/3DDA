"""
Droplet density analysis from LAMMPS dump data - Obtain ellipse of intersection z_ref plane with fitted ellipsoid
"""

# Standart python packages
import numpy as np
import math

def z_ref_ellipse(el_center, el_radii, z_ref, ntp):

    z_ref_ellipse_coeffs = get_z_ref_ellipse_coeffs(el_center, el_radii, z_ref)
    #print("z_ref_ellipse cartesian coefficients:")
    #print('a, b, c, d, e, f = ', z_ref_ellipse_coeffs)

    z_ref_ellipse_params = cart_to_pol(z_ref_ellipse_coeffs)
    #print("z_ref_ellipse parameters:")
    #print('x0, y0, ap, bp, e, phi = ', z_ref_ellipse_params)

    # Define intersection points of z_ref plane and ellipsoid: ellipse points
    x0, y0, ap, bp, e, phi = z_ref_ellipse_params
    tmin=0
    tmax=2*np.pi
    # Set a grid of the "t" parametric variable
    t = np.linspace(tmin, tmax, ntp)
    # Get ellipse points on grid t
    x_intersection_point = x0 + ap * np.cos(t) * np.cos(phi) - bp * np.sin(t) * np.sin(phi)
    y_intersection_point = y0 + ap * np.cos(t) * np.sin(phi) + bp * np.sin(t) * np.cos(phi)
    z_intersection_point = np.full_like(x_intersection_point, z_ref)
    angle_order_xy_intersection_points = np.degrees(t)
    #print("angle_order_xy_intersection_points: ", angle_order_xy_intersection_points)
    #print("Number of Intersection Points:", z_intersection_point.shape)
    #for xyz_intersection_point in zip(x_intersection_point, y_intersection_point, z_intersection_point, strict=True):
    #    print("Intersection_Point(x, y, z):", xyz_intersection_point)

    return x_intersection_point, y_intersection_point, z_intersection_point, z_ref_ellipse_params, angle_order_xy_intersection_points

def get_z_ref_ellipse_coeffs(el_center, el_radii, z_ref):
    """
    From ellipsoid equation:
       (x - el_center_x)^2/el_radii_a^2 + (y - el_center_y)^2/el_radii_b^2 + (z - el_center_z)^2/el_radii_c^2 = 1
    Get the ellipse (intersection of horizontal z_ref plane with fitted ellipsoid) equation:
       ax^2 + bxy + cy^2 + dx + ey + f = 0
    Substitute z with z_ref
    """

    a = 1/el_radii[0]**2
    b = 0
    c = 1/el_radii[1]**2
    d = -2*el_center[0]/el_radii[0]**2
    e = -2*el_center[1]/el_radii[1]**2
    f = el_center[0]**2/el_radii[0]**2 + el_center[1]**2/el_radii[1]**2 + (z_ref - el_center[2])**2/el_radii[2]**2 - 1
    
    return a, b, c, d, e, f
#
def cart_to_pol(coeffs):
    """
    Convert the cartesian conic coefficients, (a, b, c, d, e, f), to the
    ellipse parameters, where F(x, y) = ax^2 + bxy + cy^2 + dx + ey + f = 0.
    The returned parameters are x0, y0, ap, bp, e, phi, where (x0, y0) is the
    ellipse centre; (ap, bp) are the semi-major and semi-minor axes,
    respectively; e is the eccentricity; and phi is the rotation of the semi-
    major axis from the x-axis.
    """

    # We use the formulas from https://mathworld.wolfram.com/Ellipse.html
    # which assumes a cartesian form ax^2 + 2bxy + cy^2 + 2dx + 2fy + g = 0.
    # Therefore, rename and scale b, d and f appropriately.
    a = coeffs[0]
    b = coeffs[1] / 2
    c = coeffs[2]
    d = coeffs[3] / 2
    f = coeffs[4] / 2
    g = coeffs[5]

    den = b**2 - a*c
    if den > 0:
        raise ValueError('coeffs do not represent an ellipse: b^2 - 4ac must'
                         ' be negative!')

    # The location of the ellipse centre.
    x0, y0 = (c*d - b*f) / den, (a*f - b*d) / den

    num = 2 * (a*f**2 + c*d**2 + g*b**2 - 2*b*d*f - a*c*g)
    fac = np.sqrt((a - c)**2 + 4*b**2)
    # The semi-major and semi-minor axis lengths (these are not sorted).
    ap = np.sqrt(num / den / (fac - a - c))
    bp = np.sqrt(num / den / (-fac - a - c))

    # Sort the semi-major and semi-minor axis lengths but keep track of
    # the original relative magnitudes of width and height.
    width_gt_height = True
    if ap < bp:
        width_gt_height = False
        ap, bp = bp, ap

    # The eccentricity.
    r = (bp/ap)**2
    if r > 1:
        r = 1/r
    e = np.sqrt(1 - r)

    # The angle of anticlockwise rotation of the major-axis from x-axis.
    if b == 0:
        phi = 0 if a < c else np.pi/2
    else:
        phi = np.arctan((2.*b) / (a - c)) / 2
        if a > c:
            phi += np.pi/2
    if not width_gt_height:
        # Ensure that phi is the angle to rotate to the semi-major axis.
        phi += np.pi/2
    phi = phi % np.pi

    return x0, y0, ap, bp, e, phi
