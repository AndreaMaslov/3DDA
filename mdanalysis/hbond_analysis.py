"""
Droplet density analysis from LAMMPS dump data - Calculate Hydrogen bonding using MDAnalysis HydrogenBondAnalysis module
"""

# Standart python packages
import numpy as np

# MDAnalysis HydrogenBondAnalysis module
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA

def hydrogen_bonds(mda_universe, lst_number_of_hydrogen_bonds, output_csv):
    
    # Calculate Hydrogen bonding using MDAnalysis HydrogenBondAnalysis module
    dda_hba = HBA(
                universe=mda_universe,
                donors_sel='type 4',
                hydrogens_sel='type 3',
                acceptors_sel='type 1',
                d_h_cutoff=1.2,
                d_a_cutoff=2.75,
                d_h_a_angle_cutoff=140
                )
    dda_hba.run()
    dda_hba_results_hbonds = dda_hba.results.hbonds
    dda_hba_count_by_time = dda_hba.count_by_time()
    #dda_hba_count_by_ids = dda_hba.count_by_ids()
    #dda_hba_count_by_type = dda_hba.count_by_type()
    #dda_hba_tau_timeseries,  dda_hba_timeseries = dda_hba.lifetime()
    #
    # print("Statistics for dda_hba_results_hbonds ndarray:")
    # print("    Number of axes (dimensions): ", dda_hba_results_hbonds.ndim)
    # print("    Number of elements in each dimension: ", dda_hba_results_hbonds.shape)
    # print("    All dda_hba_results_hbonds elements: ", dda_hba_results_hbonds)
    #
    # print("Statistics for dda_hba_count_by_time ndarray:")
    # print("    Number of axes (dimensions): ", dda_hba_count_by_time.ndim)
    # print("    Number of elements in each dimension: ", dda_hba_count_by_time.shape)
    # print("    All dda_hba_count_by_time elements: ", dda_hba_count_by_time)
    #
    # print("Statistics for dda_hba_count_by_ids ndarray:")
    # print("    Number of axes (dimensions): ", dda_hba_count_by_ids.ndim)
    # print("    Number of elements in each dimension: ", dda_hba_count_by_ids.shape)
    # print("    All dda_hba_count_by_ids elements: ", dda_hba_count_by_ids)
    #
    # print("Statistics for dda_hba_count_by_type ndarray:")
    # print("    Number of axes (dimensions): ", dda_hba_count_by_type.ndim)
    # print("    Number of elements in each dimension: ", dda_hba_count_by_type.shape)
    # print("    All dda_hba_count_by_type elements: ", dda_hba_count_by_type)
    #
    # print("Statistics for dda_hba_tau_timeseries ndarray:")
    # print("    Number of axes (dimensions): ", dda_hba_tau_timeseries.ndim)
    # print("    Number of elements in each dimension: ", dda_hba_tau_timeseries.shape)
    # print("    All dda_hba_tau_timeseries elements: ", dda_hba_tau_timeseries)
    #

    #
    # Add to list the number of hydrogen bonds for this frame 
    lst_number_of_hydrogen_bonds.append(dda_hba_count_by_time[0])
    #

    #
    # Get distance: length of the hydrogen bond
    # Get angle: angle of the hydrogen bond
    #
    hb_distance = dda_hba_results_hbonds[:, 4]
    hb_angle = dda_hba_results_hbonds[:, 5]
    #
    # print("Statistics for hb_distance ndarray:")
    # print("    Number of axes (dimensions): ", hb_distance.ndim)
    # print("    Number of elements in each dimension: ", hb_distance.shape)
    # print("    All hb_distance elements: ", hb_distance)
    #
    # print("Statistics for hb_angle ndarray:")
    # print("    Number of axes (dimensions): ", hb_angle.ndim)
    # print("    Number of elements in each dimension: ", hb_angle.shape)
    # print("    All hb_angle elements: ", hb_angle)
    #

    #
    # Calculate the number of hydrogen bonds as a function of z position of the donor atom
    #
    hb_donor_indx = dda_hba_results_hbonds[:, 1].astype(int)  # 0-based
    hb_donor_atoms = mda_universe.atoms[hb_donor_indx]
    hb_donor_z = hb_donor_atoms.positions[:, 2]
    # print("Statistics for hb_donor_z ndarray:")
    # print("    Number of axes (dimensions): ", hb_donor_z.ndim)
    # print("    Number of elements in each dimension: ", hb_donor_z.shape)
    # print("    All hb_donor_z elements: ", hb_donor_z)
    #
    if hb_donor_z.size > 0:
        hb_donor_z_min = np.min(hb_donor_z)
        hb_donor_z_max = np.max(hb_donor_z)
        #print("hb_donor_z_min: ", hb_donor_z_min)
        #print("hb_donor_z_max: ", hb_donor_z_max)
        #
        # bins in z for the histogram
        hb_donor_z_bin_edges = np.linspace(hb_donor_z_min, hb_donor_z_max, 51)
        # print("Statistics for hb_donor_z_bin_edges ndarray:")
        # print("    Number of axes (dimensions): ", hb_donor_z_bin_edges.ndim)
        # print("    Number of elements in each dimension: ", hb_donor_z_bin_edges.shape)
        # print("    All hb_donor_z_bin_edges elements: ", hb_donor_z_bin_edges)
        #
        hb_donor_z_bin_centers = hb_donor_z_bin_edges[:-1] + 0.5
        # print("Statistics for hb_donor_z_bin_centers ndarray:")
        # print("    Number of axes (dimensions): ", hb_donor_z_bin_centers.ndim)
        # print("    Number of elements in each dimension: ", hb_donor_z_bin_centers.shape)
        # print("    All hb_donor_z_bin_centers elements: ", hb_donor_z_bin_centers)
        #
        # initialize array for counts (this is faster and more memory efficient than appending to a list)
        hb_donor_z_counts = np.full(hb_donor_z_bin_centers.size, fill_value=0.0)
        #
        for frame, donor_ix, *_ in dda_hba_results_hbonds:
            # frame variable is not used!
            donor_atom = mda_universe.atoms[donor_ix.astype(int)]
            #
            donor_zpos = donor_atom.position[2]
            hist_values, *_ = np.histogram(donor_zpos, bins = hb_donor_z_bin_edges)
            hb_donor_z_counts += hist_values * 2  # multiply by two as each hydrogen bond involves two water molecules
        #
        # print("Statistics for hb_donor_z_counts ndarray:")
        # print("    Number of axes (dimensions): ", hb_donor_z_counts.ndim)
        # print("    Number of elements in each dimension: ", hb_donor_z_counts.shape)
        # print("    All hb_donor_z_counts elements: ", hb_donor_z_counts)
        #
        # Save hbonds_donor_z.csv file
        hbonds_donor_z_ndarray = np.stack((hb_donor_z_bin_centers.round(1), hb_donor_z_counts), axis=1)
        # print("Statistics for hbonds_donor_z_ndarray:")
        # print("    Number of axes (dimensions): ", hbonds_donor_z_ndarray.ndim)
        # print("    Number of elements in each dimension: ", hbonds_donor_z_ndarray.shape)
        # print("    All hbonds_donor_z_ndarray elements: ", hbonds_donor_z_ndarray)
        np.savetxt(output_csv, hbonds_donor_z_ndarray, fmt=['%1.1f', '%u'], delimiter=';', header='z position of the donor atom; number of hydrogen bonds', comments='')
    else:
        print("For this frame there are no hydrogen bond donor atoms!")

    return lst_number_of_hydrogen_bonds
