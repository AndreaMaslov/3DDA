"""
Droplet density analysis from LAMMPS dump data - Calculate density using MDAnalysis LinearDensity module
"""

# Standart python packages
import numpy as np

# MDAnalysis LinearDensity module
from MDAnalysis.analysis.lineardensity import LinearDensity

def linear_density(group_H2O, output_csv):

    # Calculate density using MDAnalysis LinearDensity module
    dda_ldens = LinearDensity(
                group_H2O,
                binsize=0.1  # Bin width in Angstrom used to build linear density histograms.
                )
    dda_ldens.run()
    #dda_ldens_nbins = dda_ldens.nbins
    dda_ldens_results = dda_ldens.results
    #dda_ldens_results_z_keys = dda_ldens_results.z.keys()
    #print("dda_ldens_results_z_keys: ", dda_ldens_results_z_keys)
    #
    # ----- mass densities are returned in units of g/cm^3 -----
    #
    dda_ldens_results_z_mass_density = dda_ldens_results.z.mass_density
    # print("Statistics for dda_ldens_results_z_mass_density ndarray:")
    # print("    Number of axes (dimensions): ", dda_ldens_results_z_mass_density.ndim)
    # print("    Number of elements in each dimension: ", dda_ldens_results_z_mass_density.shape)
    # print("    All dda_ldens_results_z_mass_density elements: ", dda_ldens_results_z_mass_density)
    #
    dda_ldens_results_z_hist_bin_edges = dda_ldens_results.z.hist_bin_edges[:-1]
    # print("Statistics for dda_ldens_results_z_hist_bin_edges ndarray:")
    # print("    Number of axes (dimensions): ", dda_ldens_results_z_hist_bin_edges.ndim)
    # print("    Number of elements in each dimension: ", dda_ldens_results_z_hist_bin_edges.shape)
    # print("    All dda_ldens_results_z_mass_density elements: ", dda_ldens_results_z_hist_bin_edges)
    #
    # Save ldens_mass_density_z.csv file
    ldens_mass_density_z_ndarray = np.stack((dda_ldens_results_z_hist_bin_edges.round(1), dda_ldens_results_z_mass_density), axis=1)
    # print("Statistics for ldens_mass_density_z_ndarray:")
    # print("    Number of axes (dimensions): ", ldens_mass_density_z_ndarray.ndim)
    # print("    Number of elements in each dimension: ", ldens_mass_density_z_ndarray.shape)
    # print("    All ldens_mass_density_z_ndarray elements: ", ldens_mass_density_z_ndarray)
    np.savetxt(output_csv, ldens_mass_density_z_ndarray, fmt=['%1.1f', '%1.2f'], delimiter=';', header='z centers of edges of histogram bins ; mass density in the z direction', comments='')

    return
