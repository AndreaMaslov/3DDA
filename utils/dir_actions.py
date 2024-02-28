"""
Droplet density analysis from LAMMPS dump data - Directory actions module
"""

import os
import shutil

def clean(parent_dir, dir_to_clean):

    path = os.path.join(parent_dir, dir_to_clean)
    # Delete
    if os.path.exists(path):
        shutil.rmtree(path)
    # Create
    if not os.path.exists(path):
        os.mkdir(path)

def clean_new(parent_dir, dir_to_clean):

    path = os.path.join(parent_dir, dir_to_clean)
    #print("path: ", path)
    # Delete
    if os.path.exists(path):
        shutil.rmtree(path)
    # Create
    if not os.path.exists(path):
        os.makedirs(path)

def make_output_dir_from_input_dir(input_data_dir, output_data_dir):

    # In this case string BEGINS with part1 value.
    part1 = "./output/"

    # Find the index of the part1 in the string
    index = output_data_dir.index(part1)

    # Get the second part of the string
    part2 = output_data_dir[index + len(part1):]

    return part1 + input_data_dir + "/" + part2

