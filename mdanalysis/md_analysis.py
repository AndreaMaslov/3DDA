"""
Droplet density analysis from LAMMPS dump data - MDAnalysis module
"""

# Standart python packages
import numpy as np
import math

# Packages to install with "pip3 install -U MDAnalysis" (pip install -U MDAnalysisTests)
#
#   MDAnalysis
#
import MDAnalysis as mda
#
#   SciPy
#
from scipy.constants import Avogadro as av_num
from scipy.interpolate import griddata

def density(input_file_full_name, config_mda):

    #print("MDAnalysis version: ", mda.__version__)

    # Create Universe from input lammps dump file which contains ONE frame only, ergo dt is set to 1.0 ps (default value)
    mda_universe = mda.Universe(input_file_full_name, format='LAMMPSDUMP', dt=1.0)
    #
    #   Dimensions of the simulation unit cell (box) at the current time step (trajectory has ONE frame only)
    #print("mda_universe.dimensions: ", mda_universe.dimensions)
    #cell_xdim = mda_universe.dimensions[0]
    #cell_ydim = mda_universe.dimensions[1]
    #cell_zdim = mda_universe.dimensions[2]
    #print("cell_xdim, cell_ydim, cell_zdim: ", cell_xdim, cell_ydim, cell_zdim)
    #cell_diagonal = math.sqrt(cell_xdim**2 + cell_ydim**2 + cell_zdim**2)
    #print("cell_diagonal: ", cell_diagonal)
    #
    #   After Universe has been created, it has the "default AtomGroup" u.atoms which contans all loaded atoms
    #

    mda_universe.add_TopologyAttr('charge')

    # For given atom type change atom mass from guessed value 1.0 to the real value and set charge value
    #
    for atom in mda_universe.atoms:
        match atom.type:
            case '1':
                atom.mass = config_mda['type_1']['mass']
                atom.charge = config_mda['type_1']['charge']
            case '2':
                atom.mass = config_mda['type_2']['mass']
                atom.charge = config_mda['type_2']['charge']
            case '3':
                atom.mass = config_mda['type_3']['mass']
                atom.charge = config_mda['type_3']['charge']
            case '4':
                atom.mass = config_mda['type_4']['mass']
                atom.charge = config_mda['type_4']['charge']
            case '5':
                atom.mass = config_mda['type_5']['mass']
                atom.charge = config_mda['type_5']['charge']
            case '6':
                atom.mass = config_mda['type_6']['mass']
                atom.charge = config_mda['type_6']['charge']
            case _:   # Error
                print("Error: something went wrong with masses and charges!")
                atom.mass = 0.0
                atom.charge = 0.0


    # Create AtomGroup
    #
    if config_mda['OSI'] is not None:
        group_OSI = mda_universe.select_atoms(config_mda['OSI'])
    else:
        group_OSI = None
    #
    if config_mda['SI'] is not None:
        group_SI = mda_universe.select_atoms(config_mda['SI'])
    else:
        group_SI = None
    #
    if config_mda['HSur'] is not None:
        group_HSur = mda_universe.select_atoms(config_mda['HSur'])
    else:
        group_HSur = None
    #
    if config_mda['OSur'] is not None:
        group_OSur = mda_universe.select_atoms(config_mda['OSur'])
    else:
        group_OSur = None
    #
    if config_mda['H'] is not None:
        group_H = mda_universe.select_atoms(config_mda['H'])
    else:
        group_H = None
    #
    if config_mda['O'] is not None:
        group_O = mda_universe.select_atoms(config_mda['O'])
    else:
        group_O = None
    #
    if config_mda['H2O'] is not None:
        group_H2O = mda_universe.select_atoms(config_mda['H2O'])
    else:
        group_H2O = None


    #
    # Get atom positions for needed groups
    #
    sel_atoms_xyz = group_H2O.positions
    #print("Statistics for sel_atoms_xyz ndarray:")
    #print("    Number of axes (dimensions): ", sel_atoms_xyz.ndim)
    #print("    Number of elements in each dimension: ", sel_atoms_xyz.shape)
    #print("    All sel_atoms_xyz elements: ", sel_atoms_xyz)
    #
    # Get X,Y and Z coordinats from sel_atoms_xyz (slicing)
    sel_atoms_x = sel_atoms_xyz[:, 0]
    sel_atoms_y = sel_atoms_xyz[:, 1]
    sel_atoms_z = sel_atoms_xyz[:, 2]
    #
    # Get coordinates of center of mass (COM)
    sel_atoms_COM_xyz = group_H2O.center_of_mass()
    #x_com = sel_atoms_COM_xyz[0]
    #y_com = sel_atoms_COM_xyz[1]
    z_com = sel_atoms_COM_xyz[2]
    #print("Center of mass coordinates shape: ", sel_atoms_COM_xyz.shape)
    #print("Center of mass coordinates: ", sel_atoms_COM_xyz)
    #print("Center of mass coordinates (x,y,z): ", x_com, y_com, z_com)
    #
    # Calculate bounding sphere of the AtomGroup
    sel_atoms_radius_bsphere, sel_atoms_center_bsphere = group_H2O.bsphere()
    #print('sel_atoms_radius_bsphere: ', sel_atoms_radius_bsphere)
    #print('sel_atoms_center_bsphere: ', sel_atoms_center_bsphere)


    #
    # Calculate histogramdd nbins values by using MDAnalysis DensityAnalysis approach:
    #     MDAnalysis.analysis.density:
    #     https://docs.mdanalysis.org/2.4.1/_modules/MDAnalysis/analysis/density.html#DensityAnalysis
    #
    #     fixedwidth_bins:
    #     https://docs.mdanalysis.org/stable/_modules/MDAnalysis/lib/util.html#fixedwidth_bins
    #
    mda_delta = config_mda['mda_delta']    # Bin size for the density grid in ångström (same in x,y,z).
    #
    cell_vol = (mda_delta) * (mda_delta) * (mda_delta) * 1.0e-24    # Volume (in cm³) of bins voxel: cube with mda_delta edge
    #
    mda_padding = config_mda['mda_padding']   # Increase histogram dimensions by padding (on top of initial box size) in ångström.
    #
    mda_smin = np.min(sel_atoms_xyz, axis=0) - mda_padding
    mda_smax = np.max(sel_atoms_xyz, axis=0) + mda_padding
    #mda_smin = np.min(sel_atoms_xyz, axis=0)
    #mda_smax = np.max(sel_atoms_xyz, axis=0)
    #print("mda_smin, mda_smax: ", mda_smin, mda_smax)
    #print("sel_atoms_x.min(), sel_atoms_x.max():", sel_atoms_x.min(), sel_atoms_x.max())
    #print("sel_atoms_y.min(), sel_atoms_y.max():", sel_atoms_y.min(), sel_atoms_y.max())
    #print("sel_atoms_z.min(), sel_atoms_z.max():", sel_atoms_z.min(), sel_atoms_z.max())
    #
    def fixedwidth_bins(delta, xmin, xmax):
        """
        Return bins of width `delta` that cover `xmin`, `xmax` (or a larger range).
        The bin parameters are computed such that the bin size `delta` is guaranteed.
        In order to achieve this, the range `[xmin, xmax]` can be increased.
        """

        if not np.all(xmin < xmax):
            raise ValueError('Boundaries are not sane: should be xmin < xmax.')
        _delta = np.asarray(delta, dtype=np.float_)
        _xmin = np.asarray(xmin, dtype=np.float_)
        _xmax = np.asarray(xmax, dtype=np.float_)
        _length = _xmax - _xmin
        N = np.ceil(_length / _delta).astype(np.int_)  # number of bins
        dx = 0.5 * (N * _delta - _length)  # add half of the excess to each end
        
        return {'Nbins': N, 'delta': _delta, 'min': _xmin - dx, 'max': _xmax + dx}
    #
    BINS = fixedwidth_bins(mda_delta, mda_smin, mda_smax)
    mda_arange = np.transpose(np.vstack((BINS['min'], BINS['max'])))
    mda_bins = BINS['Nbins']
    #print("mda_arange: ", mda_arange)
    #print("mda_bins: ", mda_bins)
    #
    # Calculate density using numpy function histogramdd
    #
    hist, [bins_x, bins_y, bins_z] = np.histogramdd(
            (sel_atoms_x, sel_atoms_y, sel_atoms_z),
            bins = mda_bins,
            range = mda_arange,
            weights = group_H2O.masses,
            density = False
        )
    #


    #
    # Work with histogramdd output values
    #
    #print("Statistics for hist ndarray:")
    #print("    Number of axes (dimensions): ", hist.ndim)
    #print("    Number of elements in each dimension: ", hist.shape)
    #print("    Total number of elements (grid cells): ", hist.size)
    #print("    Number of non zero elements: ", np.count_nonzero(hist))
    #print("    Total number of initial data elements: ", group_H2O.n_atoms)
    #print("    Total mass of initial data elements: ", np.sum(group_H2O.masses))
    #print("    All hist elements: ", hist)
    #
    #print("sel_atoms_x.min(), sel_atoms_x.max():", sel_atoms_x.min(), sel_atoms_x.max())
    #print("sel_atoms_y.min(), sel_atoms_y.max():", sel_atoms_y.min(), sel_atoms_y.max())
    #print("sel_atoms_z.min(), sel_atoms_z.max():", sel_atoms_z.min(), sel_atoms_z.max())
    #
    #print("Statistics for edges:")
    #print("    Number of elements in each edge: ", bins_x.size, bins_y.size, bins_z.size)
    #print("    All bins_x elements: ", bins_x)
    #print("    All bins_y elements: ", bins_y)
    #print("    All bins_z elements: ", bins_z)


    #
    # Get coordinates with ONLY non-zero density values (zero values excluded)
    #
    hist_nonzero_indx, hist_nonzero_indy, hist_nonzero_indz = np.nonzero(hist)
    #print("hist_nonzero_indx: ", hist_nonzero_indx)
    #print("hist_nonzero_indy: ", hist_nonzero_indy)
    #print("hist_nonzero_indz: ", hist_nonzero_indz)
    dens_nonzero_x_full = np.array(bins_x[hist_nonzero_indx])
    dens_nonzero_y_full = np.array(bins_y[hist_nonzero_indy])
    dens_nonzero_z_full = np.array(bins_z[hist_nonzero_indz])
    cell_mass_nonzero_full = np.array(hist[hist_nonzero_indx, hist_nonzero_indy, hist_nonzero_indz])
    dens_nonzero_c_full =  cell_mass_nonzero_full / (av_num * cell_vol)  # c = color
    dens_nonzero_c_full_max = dens_nonzero_c_full.max()
    dens_nonzero_c_full_min = dens_nonzero_c_full.min()
    #print("Max Density:", dens_nonzero_c_full_max)
    #print("Min Density:", dens_nonzero_c_full_min)
    #
    #for dens_nonzero_full in zip(dens_nonzero_x_full,
    #                dens_nonzero_y_full,
    #                dens_nonzero_z_full,
    #                dens_nonzero_c_full,
    #                strict=True):
    #    print("dens_nonzero_full: x, y, z, dens_value:", dens_nonzero_full)


    #
    # Get subset of dens_nonzero_full x,y,z,c ndarrays that satisfy the filter condition
    #
    density_filter_count_min = config_mda['density_filter_count_min']
    density_filter_count_max = config_mda['density_filter_count_max']
    dens_nonzero_filtered_indx = np.where(
        (dens_nonzero_c_full >= (dens_nonzero_c_full_min * density_filter_count_min))
        &
        (dens_nonzero_c_full <= (dens_nonzero_c_full_max * density_filter_count_max))
    )
    #
    dens_nonzero_x_filtered = dens_nonzero_x_full[dens_nonzero_filtered_indx]
    dens_nonzero_y_filtered = dens_nonzero_y_full[dens_nonzero_filtered_indx]
    dens_nonzero_z_filtered = dens_nonzero_z_full[dens_nonzero_filtered_indx]
    dens_nonzero_c_filtered = dens_nonzero_c_full[dens_nonzero_filtered_indx]
    #
    #for dens_nonzero_filtered in zip(dens_nonzero_x_filtered,
    #                dens_nonzero_y_filtered,
    #                dens_nonzero_z_filtered,
    #                dens_nonzero_c_filtered,
    #                strict=True):
    #    print("dens_nonzero_filtered: x, y, z, dens:", dens_nonzero_filtered)


    #
    # Create a dense multi-dimensional "meshgrid" array from a set of input arrays
    #
    #   Density points coordinates
    dens_points = np.array(
        (
            dens_nonzero_x_filtered.flatten(),
            dens_nonzero_y_filtered.flatten(),
            dens_nonzero_z_filtered.flatten()
        )
    ).T
    #
    #   Density values @ above density points
    dens_values = dens_nonzero_c_filtered.flatten()
    #
    #   Create mesh grid
    mgrid_x, mgrid_y, mgrid_z = np.mgrid[
        dens_nonzero_x_filtered.min():dens_nonzero_x_filtered.max():(bins_x.size * 1j),
        dens_nonzero_y_filtered.min():dens_nonzero_y_filtered.max():(bins_y.size * 1j),
        dens_nonzero_z_filtered.min():dens_nonzero_z_filtered.max():(bins_z.size * 1j)
        ]
    mgrid_xyz = np.column_stack([mgrid_x.flatten(), mgrid_y.flatten(), mgrid_z.flatten()])
    #
    #   Interpolate density values from dens_points to nodes of mesh grid
    mgrid_density = griddata(dens_points, dens_values, (mgrid_x ,mgrid_y, mgrid_z), method='linear', fill_value=0.0)
    #
    #print("Statistics for mgrid_x:")
    #print("    Number of axes (dimensions): ", mgrid_x.ndim)
    #print("    Number of elements in each dimension: ", mgrid_x.shape)
    #print("    Total number of elements: ", mgrid_x.size)
    #print("    All mgrid_x elements: ", mgrid_x)
    #
    # print("Statistics for mgrid_x.flatten():")
    # print("    Number of axes (dimensions): ", mgrid_x.flatten().ndim)
    # print("    Number of elements in each dimension: ", mgrid_x.flatten().shape)
    # print("    Total number of elements: ", mgrid_x.flatten().size)
    #print("    All mgrid_x.flatten() elements: ", mgrid_x.flatten())
    #
    # print("Statistics for mgrid_density:")
    # print("    Number of axes (dimensions): ", mgrid_density.ndim)
    # print("    Number of elements in each dimension: ", mgrid_density.shape)
    # print("    Total number of elements: ", mgrid_density.size)
    #print("    All mgrid_density elements: ", mgrid_density)
    #
    # print("Statistics for mgrid_density.flatten():")
    # print("    Number of axes (dimensions): ", mgrid_density.flatten().ndim)
    # print("    Number of elements in each dimension: ",mgrid_density.flatten().shape)
    # print("    Total number of elements: ", mgrid_density.flatten().size)
    #print("    All mgrid_density.flatten() elements: ", mgrid_density.flatten())


    #
    # Prepare for Convex hull of xyz nonzero density points
    #
    #   Only nonzero density point give the meaningful convex hull form.
    #   Otherwise convex hull will be always parallelepiped, because all external point have zero density!
    #
    dens_nonzero_xyz_filtered = np.hstack((dens_nonzero_x_filtered[:,np.newaxis], dens_nonzero_y_filtered[:,np.newaxis], dens_nonzero_z_filtered[:,np.newaxis]))
    #print("Statistics for dens_nonzero_xyz_filtered ndarray:")
    #print("    Number of axes (dimensions): ", dens_nonzero_xyz_filtered.ndim)
    #print("    Number of elements in each dimension: ", dens_nonzero_xyz_filtered.shape)
    #print("    All dens_nonzero_xyz_filtered elements: ", dens_nonzero_xyz_filtered)


    #
    # Select the points in the grid that have non-zero data values  
    #nonzero_dens_points = dens_points[mgrid_density != 0]

    return mda_universe, z_com, mgrid_x, mgrid_y, mgrid_z, mgrid_density, dens_nonzero_xyz_filtered, dens_nonzero_x_filtered, dens_nonzero_y_filtered, dens_nonzero_z_filtered, dens_nonzero_c_filtered, sel_atoms_radius_bsphere, sel_atoms_center_bsphere, group_H2O, group_O, group_H, group_OSI, group_SI, group_HSur, group_OSur
