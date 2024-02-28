"""
Droplet density analysis from LAMMPS dump data - Most probable angle per frame module
"""

# Standart python packages
import numpy as np

# Packages to install with "pip3 install -U package-name"
#
#   KDE
#
from sklearn.neighbors import KernelDensity
from sklearn.model_selection import GridSearchCV

def kde(lst_frames, frame_number, lst_tangent_angles_per_frame, lst_angles_mprob_values):

    # Calculate KDE of lst_angles_per_frame values with sklearn
    #   In sklearn terms lst_angles_per_frame have single feature (angle) and many (len(lst_angles_per_frame)) samples (angles's values).
    #   Fitting a model requires a 2D array. i.e (n_samples, n_features).
    #   That's why, to fit lst_angles_per_frame to KernelDensity sklearn model it should be converted to ndarray and reshaped:
    angles_per_frame = np.array(lst_tangent_angles_per_frame).reshape(-1, 1)
    #print("Statistics for angles_per_frame ndarray:")
    #print("    Number of axes (dimensions): ", angles_per_frame.ndim)
    #print("    Number of elements in each dimension: ", angles_per_frame.shape)
    #print("    Min and Max values of angles_per_frame: ", angles_per_frame.min(), angles_per_frame.max())
    #number_of_samples = angles_per_frame.shape[0]
    #print(f"Frame {frame_number} has number of samples: {number_of_samples}")

    # Use grid search cross-validation to optimize the bandwidth and get the best kernel estimator
    params = {"bandwidth": np.logspace(-1, 1, 100)}
    grid = GridSearchCV(KernelDensity(), params,)
    grid.fit(angles_per_frame)
    best_kde = grid.best_estimator_
    #print("Best Kernel: ", best_kde.kernel)
    #print("Best Bandwidth: ", best_kde.bandwidth)
    #print('accuracy =', grid.best_score_)

    # Use the best estimator to compute the kernel density estimate
    #x_samples = np.linspace(angles_per_frame.min(), angles_per_frame.max(), 5000)[:, np.newaxis]
    x_samples = np.linspace(0, 180, 50000)[:, np.newaxis]
    #print("Statistics for x_samples ndarray:")
    #print("    Number of axes (dimensions): ", x_samples.ndim)
    #print("    Number of elements in each dimension: ", x_samples.shape)
    #print("    All x_samples elements: ", x_samples)
    angles_per_frame_dens = np.exp(best_kde.score_samples(x_samples))
    #print("Statistics for angles_per_frame_dens ndarray:")
    #print("    Number of axes (dimensions): ", angles_per_frame_dens.ndim)
    #print("    Number of elements in each dimension: ", angles_per_frame_dens.shape)
    #print("    All ref_tri_angles_dens elements: ", angles_per_frame_dens)

    # Get the most probable angle value
    #print("Statistics for x_samples[:,0] ndarray:")
    #print("    Number of axes (dimensions): ", x_samples[:,0].ndim)
    #print("    Number of elements in each dimension: ", x_samples[:,0].shape)
    #print("    All x_samples[:,0] elements: ", x_samples[:,0])
    angles_per_frame_dens_max_value_indx = np.argmax(angles_per_frame_dens)
    angle_per_frame_mprob_value = x_samples[:,0][angles_per_frame_dens_max_value_indx]
    #
    print(frame_number + " ", angle_per_frame_mprob_value)
    #
    lst_frames.append(frame_number)
    lst_angles_mprob_values.append(angle_per_frame_mprob_value)
    #

    return lst_frames, lst_angles_mprob_values, x_samples, angles_per_frame_dens, angle_per_frame_mprob_value
