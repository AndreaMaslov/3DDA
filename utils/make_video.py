"""
Droplet density analysis from LAMMPS dump data - Make video module
"""

# Standart python packages
import os

# Packages to install with "pip3 install -U package-name"
#
#   Imageio for video creation
#
import imageio.v2 as iio

def make_video(images_for_video, output_gif):

    video_writer = iio.get_writer(output_gif, mode='I', duration=0.01, loop=1)

    # Iterate video images directory (presume it contains only files!)
    with os.scandir(images_for_video) as direntry:
        direntries = list(direntry)
    # Sort direntries by file name
    direntries.sort(key=lambda x: int(os.path.splitext(x.name)[0]))

    for image_file in direntries:
        image_file_full_name = os.path.abspath(image_file)
        #print('image_file_full_name: ', image_file_full_name)
        video_writer.append_data(iio.imread(image_file_full_name))
    
    video_writer.close()
