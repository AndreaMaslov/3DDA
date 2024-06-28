# 3DDA
3DDA: A novel Python toolkit to analyze 3D-Dynamic Contact Angles from Molecular Dynamics simulations

This Python script is designed to analyze the 3d contact angles of droplets from MD simulations.  In this version the custom input strucutre comes from LAMMPS trajectories, but by changing the parser it is possible to analyze trajectories from others common MD solvers.


## Usage

To run the software, use the following command:

python 3DDA.py [options]

For detailed information about the options, you can use the help flag:

python3 3DDA.py -h

Options
The script accepts the following options:

-h, --help: Show this help message and exit.
-pff_i int: Input frames processing frequency (default: 1).
-pff_o int: Output frames processing frequency (default: 1).
-zv int: Z up move value (default: 5).
-chi int: Convex hull iterations (default: 2).
-ntp int: Number of tangent points for angles calculations (default: 180)

**Requirements**
Ensure you have the following Python packages installed:

numpy
scipy
scikit learn
Plotly

**Installation**
Clone the repository and navigate to the directory:

git clone https://github.com/AndreaMaslov/3DDA
cd 3DDA

**Input Files**
Place all necessary input files in the input folder within the project directory. The script will process these files based on the provided options.

**Example**
To run the script with default options, use:
python3 3DDA.py

To run the script with custom parameters, for example:

python3 3DDA.py -pff_i 2 -pff_o 2 -zv 10 -chi 3 -ntp 200
