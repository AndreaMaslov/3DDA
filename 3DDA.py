"""
This is the entry point (main module) which imports all other modules and should be run as:
    python droplet_lammps_density_3D.py [-h] [-pff_i int] [-pff_o int] [-zv int] [-chi int] [-ntp int]
"""

# Import standard python packages
#    Python has separate concepts of importing and loading a module:
#      - Loading actually creates an in-memory representation of the module and produces a module object.
#      - Importing only binds a module in the current namespace and may load the module if it is not present
#    Example:
#      When numpy is imported for the first time, it is loaded.
#      When numpy is imported for the next time in any child module, it is binded to the already loaded numpy module.
#
import argparse
import os
import sys as syst
import numpy as np
#from time import perf_counter

# Packages to install with "pip3 install -U package-name"
#
#   PyYAML
#
import yaml

# Import local modules
import utils.dir_actions as dirs
#import utils.make_video as video
#
import plots.color_density as plt_cd
import plots.ch_fitted_ellipsoid as plt_cfe
import plots.angles_dens as plt_adn
import plots.polar_chart_per_frame as plt_pcf
import plots.tangent_angle_per_frame as plt_taf
import plots.radius_per_frame as plt_raf
import plots.area_per_frame as plt_arf
import plots.volume_per_frame as plt_vof
#
import calculation.coordinate_systems as cs
import calculation.ellipsoid_fit_convex_hull as efit
import calculation.z_ref_ellipsiod_intersection_points as zre
import calculation.z_ref_ellipsiod_tangent_angle as zeta
import calculation.most_probable_angle_per_frame as mpaf
import calculation.droplet_additional_values_per_frame as davf
#
import mdanalysis.md_analysis as mda
import mdanalysis.hbond_analysis as mda_hbond
import mdanalysis.ld_analysis as mda_ld
#
import convex_hull.density_convex_hull as dch

# Get input parameters
parser = argparse.ArgumentParser(description="Droplet density analysis from LAMMPS dump data.",
                                formatter_class=argparse.MetavarTypeHelpFormatter)
parser.add_argument("-pff_i", dest='pff_i', default=1, type=int,
                    help="input frames processing frequency (default: 1)")
parser.add_argument("-pff_o", dest='pff_o', default=1, type=int,
                    help="output frames processing frequency (default: 1)")
parser.add_argument("-zv", dest='zv', default=5, type=int,
                    help="z up move value (default: 5)")
parser.add_argument("-chi", dest='chi', default=2, type=int,
                    help="convex hull iterations (default: 2)")
parser.add_argument("-ntp", dest='ntp', default=360, type=int,
                    help="number of tangent points for angles calculations (default: 100)")
args = parser.parse_args()
print(f"Input parametrs values - pff_i: {args.pff_i}, pff_o: {args.pff_o}, zv: {args.zv}, chi: {args.chi}, ntp: {args.ntp}")
# Get config parameters
#
#   Config main module (DDA)
config_file_dda = "./config/config_dda.yaml"
with open(config_file_dda, 'r',  encoding="utf-8") as file:
    config_dda = yaml.safe_load(file)
#
input_data_dir = config_dda['input_dir'].split("/")[-1]  #  last element has an index of -1
#
output_dir_mass_density = dirs.make_output_dir_from_input_dir(input_data_dir, config_dda['output_dir_mass_density'])
output_dir_convex_hull_fit_ellipsoid = dirs.make_output_dir_from_input_dir(input_data_dir, config_dda['output_dir_convex_hull_fit_ellipsoid'])
output_dir_angles_dens = dirs.make_output_dir_from_input_dir(input_data_dir, config_dda['output_dir_angles_dens'])
output_dir_polar_chart = dirs.make_output_dir_from_input_dir(input_data_dir, config_dda['output_dir_polar_chart'])
output_dir_droplet_values_per_frame = dirs.make_output_dir_from_input_dir(input_data_dir, config_dda['output_dir_droplet_values_per_frame'])
output_dir_csv_files_all_frames = dirs.make_output_dir_from_input_dir(input_data_dir, config_dda['output_dir_csv_files_all_frames'])
output_dir_csv_files_single_frame = dirs.make_output_dir_from_input_dir(input_data_dir, config_dda['output_dir_csv_files_single_frame'])
#
#   Config MDAnalysis
#config_file_mda = "./config/config_mda_V1.yaml"
config_file_mda = "./config/config_mda_V2.yaml"
with open(config_file_mda, 'r',  encoding="utf-8") as file:
    config_mda = yaml.safe_load(file)

# Start the time counter
#time_start = perf_counter()

# Set no truncation mode for printing long arrays
# In this case use command: python droplet_lammps_density.py > test
np.set_printoptions(threshold=syst.maxsize)   # sys module should be imported

# Clean output directories
dirs.clean_new(output_dir_mass_density, config_dda['dir_html'])
dirs.clean_new(output_dir_mass_density, config_dda['dir_images'])
dirs.clean_new(output_dir_mass_density, config_dda['dir_gif'])
dirs.clean_new(output_dir_convex_hull_fit_ellipsoid, config_dda['dir_html'])
dirs.clean_new(output_dir_convex_hull_fit_ellipsoid, config_dda['dir_images'])
dirs.clean_new(output_dir_angles_dens, config_dda['dir_html'])
dirs.clean_new(output_dir_angles_dens, config_dda['dir_images'])
dirs.clean_new(output_dir_polar_chart, config_dda['dir_html'])
dirs.clean_new(output_dir_polar_chart, config_dda['dir_images'])
dirs.clean_new(output_dir_droplet_values_per_frame, "")
dirs.clean_new(output_dir_csv_files_all_frames, "")
dirs.clean_new(output_dir_csv_files_single_frame, "")
dirs.clean_new(output_dir_csv_files_single_frame, config_dda['dir_csv_hbonds'])
dirs.clean_new(output_dir_csv_files_single_frame, config_dda['dir_csv_ldens'])

# Iterate input directory with lammps dump data files (presume it contains only files!)
with os.scandir(config_dda['input_dir']) as direntry:
    direntries = list(direntry)
# Sort direntries by file name (direntries is list!)
direntries.sort(key=lambda x: int(x.name))
#print("direntries: ", direntries)
#
# Use the extended slicing syntax list[start:stop:step] to get a new list containing every nth element.
# Leave start and stop empty, and set step to the desired n.
#   Set input data list to process
input_files = direntries[::args.pff_i]
#print("input_files: ", input_files)
output_plots = input_files[::args.pff_o]
#print("output_plots: ", output_plots)

# Allocate lists for frames and corresponding most probable values of angles between ellipsoid tangent plane and reference horizontal plane in intersection points
lst_frames = []
lst_angles_mprob_values = []
lst_radius_values = []
lst_volume_values = []
lst_area_values = []
lst_number_of_hydrogen_bonds = []

# Process sorted data files using list extended slicing
for data_file in input_files:
    input_file_full_name = os.path.abspath(data_file)
    frame_number = data_file.name
    #print("data_file_full_name: ", input_file_full_name)
    #print("frame_number: ", frame_number)

    # Use MDAnalysis to obtain droplet non-zero density points and relative COM (center of masses)
    mda_universe, z_com, mgrid_x, mgrid_y, mgrid_z, mgrid_density, dens_nonzero_xyz_full, dens_nonzero_x, dens_nonzero_y, dens_nonzero_z, dens_nonzero_c, sel_atoms_radius_bsphere, sel_atoms_center_bsphere, group_H2O, group_O, group_H, group_OSI, group_SI, group_HSur, group_OSur = mda.density(input_file_full_name, config_mda)
    #print("z_com: ", z_com)

    # Use HydrogenBondAnalysis for hydrogen bonding calculation
    output_csv = os.path.join(output_dir_csv_files_single_frame, config_dda['dir_csv_hbonds'], frame_number + ".csv")
    lst_number_of_hydrogen_bonds = mda_hbond.hydrogen_bonds(mda_universe, lst_number_of_hydrogen_bonds, output_csv)

    # Use LinearDensity for hydrogen bonding calculation
    output_csv = os.path.join(output_dir_csv_files_single_frame, config_dda['dir_csv_ldens'], frame_number + ".csv")
    mda_ld.linear_density(group_H2O, output_csv)

    # Obtain reference horizontal plane
    z_ref = cs.z_ref(dens_nonzero_xyz_full, args.zv)
        
    # Perform Convex Hull of 3D nonzero density points
    dens_nonzero_xyz = dens_nonzero_xyz_full
    for i in range(args.chi):
        #print("convex hull interation: ", i)
        hull_vertices_indx, hull_vertices_xyz = dch.density_convex_hull(dens_nonzero_xyz, z_ref, z_com)
        dens_nonzero_xyz = np.delete(dens_nonzero_xyz, hull_vertices_indx, axis=0)
    #print("Statistics for hull_vertices_xyz ndarray:")
    #print("    Number of axes (dimensions): ", hull_vertices_xyz.ndim)
    #print("    Number of elements in each dimension: ", hull_vertices_xyz.shape)
    #print("    All hull_vertices_xyz elements: ", hull_vertices_xyz)

    # Fit ellipsoid to Convex Hull points:
    #    (x - el_center_x)^2/el_radii_a^2 + (y - el_center_y)^2/el_radii_b^2 + (z - el_center_z)^2/el_radii_c^2 = 1
    el_center, el_radii = efit.ellipsoid_fit(hull_vertices_xyz, sel_atoms_radius_bsphere, sel_atoms_center_bsphere)
    #el_center, el_radii = efit.ellipsoid_fit(dens_nonzero_xyz_full, sel_atoms_radius_bsphere, sel_atoms_center_bsphere)

    # Obtain intersection points of z_ref plane with fitted ellipsoid. These points form ellipse
    x_intersection_point, y_intersection_point, z_intersection_point, z_ref_ellipse_params, angle_order_xy_intersection_points = zre.z_ref_ellipse(el_center, el_radii, z_ref, args.ntp)

    # Calculate angels between ellipsoid tangent plane and reference horizontal plane in intersection points
    lst_tangent_angles_per_frame = zeta.ellipsoid_tangent_angle(frame_number, el_center, el_radii, x_intersection_point, y_intersection_point, z_intersection_point)

    # Calculate KDE to obtained the most probable angle value
    lst_frames, lst_angles_mprob_values, frame_x_samples, frame_angles_dens, angle_per_frame_mprob_value= mpaf.kde(lst_frames, frame_number, lst_tangent_angles_per_frame, lst_angles_mprob_values)

    # Calculate Droplet Additional Values
    lst_radius_values, lst_area_values, lst_volume_values = davf.droplet_additional_values(frame_number, z_ref_ellipse_params, el_center, el_radii, z_ref, lst_radius_values, lst_area_values, lst_volume_values)

    # Plot results for current frame
    if output_plots.count(data_file) > 0:
        #
        output_html = os.path.join(output_dir_mass_density, config_dda['dir_html'], frame_number + ".html")
        output_png = os.path.join(output_dir_mass_density, config_dda['dir_images'], frame_number + ".png")
        plt_cd.color_density(output_html, output_png, frame_number, mgrid_x, mgrid_y, mgrid_z, mgrid_density)
        #
        output_html = os.path.join(output_dir_convex_hull_fit_ellipsoid, config_dda['dir_html'], frame_number + ".html")
        output_png = os.path.join(output_dir_convex_hull_fit_ellipsoid, config_dda['dir_images'], frame_number + ".png")
        plt_cfe.fitted_ellipsoid(output_html, output_png, frame_number, dens_nonzero_xyz_full, hull_vertices_xyz, z_ref, z_ref_ellipse_params, el_center, el_radii, x_intersection_point, y_intersection_point, z_intersection_point)
        #
        output_html = os.path.join(output_dir_angles_dens, config_dda['dir_html'], frame_number + ".html")
        output_png = os.path.join(output_dir_angles_dens, config_dda['dir_images'], frame_number + ".png")
        plt_adn.angles_dens(output_html, output_png, frame_number, frame_x_samples, frame_angles_dens)
        #
        output_html = os.path.join(output_dir_polar_chart, config_dda['dir_html'], frame_number + ".html")
        output_png = os.path.join(output_dir_polar_chart, config_dda['dir_images'], frame_number + ".png")
        plt_pcf.polar_chart(output_html, output_png, frame_number, angle_per_frame_mprob_value, lst_tangent_angles_per_frame, angle_order_xy_intersection_points)


# Display original lists
# print("lst_frames:\n", lst_frames, "\n")


# Convert lists to numpy arrays
frames_values = np.array(list(map(int, lst_frames)))
hydrogen_bonds_values = np.array(list(map(int, lst_number_of_hydrogen_bonds)))
angles_mprob_values = np.array(list(map(round, lst_angles_mprob_values)))
radius_values = np.array(list(map(round, lst_radius_values)))
area_values = np.array(list(map(round, lst_area_values)))
volume_values = np.array(list(map(round, lst_volume_values)))


# Save frames_hydrogen_bonds.csv file
frames_hydrogen_bonds_ndarray = np.stack((frames_values, hydrogen_bonds_values), axis=1)
# print("Statistics for frames_hydrogen_bonds_ndarray:")
# print("    Number of axes (dimensions): ", frames_hydrogen_bonds_ndarray.ndim)
# print("    Number of elements in each dimension: ", frames_hydrogen_bonds_ndarray.shape)
# print("    All frames_hydrogen_bonds_ndarray elements: ", frames_hydrogen_bonds_ndarray)
frames_hydrogen_bonds_csv_file = os.path.join(output_dir_csv_files_all_frames, "frames_hydrogen_bonds.csv")
np.savetxt(frames_hydrogen_bonds_csv_file, frames_hydrogen_bonds_ndarray, fmt=['%s', '%s'], delimiter=';', header='frame; number of hydrogen bonds', comments='')


# Save frames_angles.csv file
frames_angles_ndarray = np.stack((frames_values, angles_mprob_values), axis=1)
# print("Statistics for frames_angles_ndarray:")
# print("    Number of axes (dimensions): ", frames_angles_ndarray.ndim)
# print("    Number of elements in each dimension: ", frames_angles_ndarray.shape)
# print("    All frames_angles_ndarray elements: ", frames_angles_ndarray)
frames_angles_csv_file = os.path.join(output_dir_csv_files_all_frames, "frames_angles.csv")
np.savetxt(frames_angles_csv_file, frames_angles_ndarray, fmt=['%s', '%s'], delimiter=';', header='frame; angle', comments='')


# Draw in 2D the most probable angle value for every frame
output_html = os.path.join(output_dir_droplet_values_per_frame, "droplet_tangent_angle_per_frame_lines.html")
output_png = os.path.join(output_dir_droplet_values_per_frame, "droplet_tangent_angle_per_frame_lines.png")
plt_taf.tangent_angle_value_per_frame(frames_values, angles_mprob_values, output_html, output_png, 'lines')
#
output_html = os.path.join(output_dir_droplet_values_per_frame, "droplet_tangent_angle_per_frame_markers.html")
output_png = os.path.join(output_dir_droplet_values_per_frame, "droplet_tangent_angle_per_frame_markers.png")
plt_taf.tangent_angle_value_per_frame(frames_values, angles_mprob_values, output_html, output_png, 'markers')
#
# Draw in 2D the radius value for every frame
output_html = os.path.join(output_dir_droplet_values_per_frame, "droplet_radius_per_frame.html")
output_png = os.path.join(output_dir_droplet_values_per_frame, "droplet_radius_per_frame.png")
plt_raf.radius_value_per_frame(frames_values, radius_values, output_html, output_png)
#
# Draw in 2D the area value for every frame
output_html = os.path.join(output_dir_droplet_values_per_frame, "droplet_area_per_frame.html")
output_png = os.path.join(output_dir_droplet_values_per_frame, "droplet_area_per_frame.png")
plt_arf.area_value_per_frame(frames_values, area_values, output_html, output_png)
#
# Draw in 2D the volume value for every frame
output_html = os.path.join(output_dir_droplet_values_per_frame, "droplet_volume_per_frame.html")
output_png = os.path.join(output_dir_droplet_values_per_frame, "droplet_volume_per_frame.png")
plt_vof.volume_value_per_frame(frames_values, volume_values, output_html, output_png)
#


# Create video files
images_for_video = os.path.join(output_dir_mass_density, config_dda['dir_images'])
output_gif = os.path.join(output_dir_mass_density, config_dda['dir_gif'], "mass_density.gif")
#video.make_video(images_for_video, output_gif)


# Stop the time counter
#time_stop = perf_counter()
#print("Execution elapsed time in minutes: ", (time_stop-time_start)/60)
