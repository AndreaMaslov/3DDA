"""
Droplet density analysis from LAMMPS dump data - Fit ellipsoid
"""

# Standart python packages
import numpy as np
import math

# Packages to install with pip install -U "jax[cpu]"
#
#   JAX
#
import jax.numpy as jnp
from jax import grad
#from jax import random
from jax.config import config
#
#   Optimization
#
#from scipy.optimize import fmin_bfgs
from scipy.optimize import fmin_l_bfgs_b


def ellipsoid_fit(hull_vertices_xyz, sel_atoms_radius_bsphere, sel_atoms_center_bsphere):
    """
    Based on: https://jekel.me/2021/A-better-way-to-fit-Ellipsoids/

    """
    config.update('jax_enable_x64', True)
    config.update('jax_platform_name', 'cpu')
    #key = random.PRNGKey(0)

    x = hull_vertices_xyz[:, 0]
    y = hull_vertices_xyz[:, 1]
    z = hull_vertices_xyz[:, 2]

    gamma_guess = np.random.random(6)
    gamma_guess[0] = sel_atoms_center_bsphere[0]
    gamma_guess[1] = sel_atoms_center_bsphere[1]
    gamma_guess[2] = sel_atoms_center_bsphere[2]
    gamma_guess[3] = sel_atoms_radius_bsphere
    gamma_guess[4] = sel_atoms_radius_bsphere
    gamma_guess[5] = sel_atoms_radius_bsphere
    #
    #print("Initial values (from bsphere): ", gamma_guess)

    def predict(gamma):
        # compute f hat
        x0 = gamma[0]
        y0 = gamma[1]
        z0 = gamma[2]
        a2 = gamma[3]**2
        b2 = gamma[4]**2
        c2 = gamma[5]**2
        zeta0 = (x - x0)**2 / a2
        zeta1 = (y - y0)**2 / b2
        zeta2 = (z - z0)**2 / c2
        return zeta0 + zeta1 + zeta2
    #
    def loss(g):
        # compute mean squared error
        pred = predict(g)
        target = jnp.ones_like(pred)
        mse = jnp.square(pred-target).mean()
        return mse
    #
    #initial_mean_squared_error = loss(gamma_guess)
    #print("initial_mean_squared_error: ", initial_mean_squared_error)

    # Compute the derivatives of the mean squared error with respect to each gamma component
    #deriv_mse_gamma_components = grad(loss)(gamma_guess)
    #print("deriv_mse_gamma_components: ", deriv_mse_gamma_components)

    # Optimize (minimize) error function
    el_bounds = [(None, None),
                 (None, None),
                 (None, None),
                 (- sel_atoms_radius_bsphere, None),
                 (- sel_atoms_radius_bsphere, None),
                 (- sel_atoms_radius_bsphere, + sel_atoms_radius_bsphere)]
    minimized_function = fmin_l_bfgs_b(
        loss,
        gamma_guess,
        args=(),
        bounds=el_bounds,
        fprime=grad(loss),
        disp=0,  # disp=1,
        callback=None
    )

    #print("minimized_function: ", minimized_function)
    #final_mean_squared_error = minimized_function[1]
    #print("Fitting error: ", final_mean_squared_error)
    ellipsoid_parameters = minimized_function[0]
    #print("ellipsoid_parameters: ", ellipsoid_parameters)

    center = [ellipsoid_parameters[0], ellipsoid_parameters[1], ellipsoid_parameters[2]] # Ellipsoid center(x,y,z)
    radii = [ellipsoid_parameters[3], ellipsoid_parameters[4], ellipsoid_parameters[5]]  # Ellipsoid radii(a,b,c)

    return center, radii
