import numpy as np
import fractions
import skimage.transform as skit
from typing import Tuple
from ct_tools.ctdata import CTdata


def resample_ctdata(ct: CTdata, voxels: Tuple[float, float, float], size_or_voxels='size') -> CTdata:
    # based on the the resampleCT subroutine in ctcreate.mortran in EGSnrc
    # resampling weighs the original CT values by the fractional volume of each resampled (xyz) voxel that overlaps
    # with a given CT voxel
    if size_or_voxels == "size":
        voxel_size_in_cm = voxels
    elif size_or_voxels == "voxels":
        voxel_size_in_cm = [ct.image_size_in_cm()[i] / v for (i, v) in enumerate(voxels)]
    else:
        raise Exception("Argument 'size_or_count' must be either 'size' or 'voxels'.")
    print("Resampling CT data...")
    print("Requested voxel size:\n{:.4f} mm x {:.4f} mm x {:.4f} mm\n".format(voxel_size_in_cm[0] * 10,
                                                                              voxel_size_in_cm[1] * 10,
                                                                              voxel_size_in_cm[2] * 10))

    check_requested_voxel_size(ct.image_size_in_cm(), voxel_size_in_cm)
    dimensions, adjusted_voxel_size = adjust_requested_voxel_size(ct.image_size_in_cm(), voxel_size_in_cm)
    bounds = calculate_new_bounds(ct.bounds, adjusted_voxel_size, dimensions)

    resampled = CTdata()
    resampled.dimensions = dimensions
    resampled.bounds = bounds
    resampled.image = np.zeros(dimensions)

    resampled.image = skit.resize(ct.image, dimensions, preserve_range=True)

    return resampled


def check_requested_voxel_size(ct_image_size, voxel_size):
    if any(v < 0 for v in voxel_size):
        raise Exception("Voxel size must be positive.")

    if any(v > ct_image_size[i] for (i, v) in enumerate(voxel_size)):
        raise Exception("Voxel size in any one direction cannot be greater"
                        "than the size of the original CT in this direction.")


def adjust_requested_voxel_size(ct_image_size, voxel_size):
    dimensions = [int(ct_image_size[i] / v) for (i, v) in enumerate(voxel_size)]
    adjusted_voxel_size = [float(ct_image_size[i] / n) for (i, n) in enumerate(dimensions)]

    print("Adjusted dimensions so that an integer number of voxels fit exactly on the CT data:\n")
    print("New voxel size:\n{:.4f} mm x {:.4f} mm x {:.4f} mm\n".format(adjusted_voxel_size[0] * 10,
                                                                        adjusted_voxel_size[1] * 10,
                                                                        adjusted_voxel_size[2] * 10))
    print("New number of voxels:\n{:d} x {:d} x {:d}\n".format(dimensions[0],
                                                               dimensions[1],
                                                               dimensions[2]))
    return dimensions, adjusted_voxel_size


def calculate_new_bounds(ct_bounds, voxel_size, dimensions):
    print("Calculating new bounds...\n")
    return [[ct_bounds[i][0] + j * voxel_size[i] for j in range(dimensions[i] + 1)] for i in range(3)]


def fast_downscale(old_size, new_size):
    original_dim = old_size
    new_dim = new_size
    f = fractions.Fraction(original_dim, new_dim).limit_denominator()
    new_dim_reduced = f.denominator
    original_dim_reduced = f.numerator
    weights = np.zeros((new_dim_reduced, original_dim_reduced))
    wtot = original_dim / new_dim
    index_orig = 0
    index = 0
    wrem = wtot

    while index < new_dim_reduced and index_orig < original_dim_reduced:
        while wrem >= 1.:
            weights[(index, index_orig)] = 1
            index_orig += 1
            wrem -= 1
        if index_orig == original_dim_reduced:
            break
        weights[(index, index_orig)] = wrem
        next_weight = 1 - wrem
        index += 1
        if index == new_dim_reduced:
            break
        weights[(index, index_orig)] = next_weight
        wrem = wtot - next_weight
        index_orig += 1

    the_weights = np.zeros((new_dim, original_dim))
    n = int(original_dim / original_dim_reduced)
    for i in range(n):
        np.copyto(the_weights[i*new_dim_reduced:(i+1)*new_dim_reduced,
                  i*original_dim_reduced:(i+1)*original_dim_reduced], weights)

    return the_weights / wtot

