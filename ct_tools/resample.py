import numpy as np
from typing import Tuple
from ct_tools.ctdata import CTdata


def resample_ctdata(ct: CTdata, voxel_size_in_cm: Tuple[float, float, float]) -> CTdata:
    resampled = CTdata()

    check_requested_voxel_size(ct.image_size_in_cm(), voxel_size_in_cm)

    dimensions, adjusted_voxel_size = adjust_requested_voxel_size(ct.image_size_in_cm(), voxel_size_in_cm)

    bounds = calculate_new_bounds(ct.bounds, adjusted_voxel_size, dimensions)

    ct_voxel_size = ct.voxel_size_in_cm()
    for (ct_ijk, value) in np.ndenumerate(ct.image):
        lower_ijk, upper_ijk = find_indices_of_resampled_voxels_where_boundaries_of_ct_voxel_lie(ct.bounds, bounds,
                                                                                                 ct_ijk)


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

    print("Adjusted dimensions (in cm) so that an integer number of voxels fit exactly on the CT data:")
    print("New voxel size in cm: ({:.2f}, {:.2f}, {:.2f})".format(adjusted_voxel_size[0],
                                                                  adjusted_voxel_size[1],
                                                                  adjusted_voxel_size[2]))
    print("New number of voxels: ({:d}, {:d}, {:d})".format(dimensions[0],
                                                            dimensions[1],
                                                            dimensions[2]))
    return dimensions, adjusted_voxel_size


def calculate_new_bounds(ct_bounds, voxel_size, dimensions):
    return [[ct_bounds[i][0] + j * voxel_size[i] for j in range(dimensions[i] + 1)] for i in range(3)]


def find_indices_of_resampled_voxels_where_boundaries_of_ct_voxel_lie(ct_bounds, bounds, ct_ijk):
    ct_lower_bounds = [ct_bounds[i][v] for (i, v) in enumerate(ct_ijk)]
    ct_upper_bounds = [ct_bounds[i][v + 1] for (i, v) in enumerate(ct_ijk)]

    lower_indices = np.array([np.searchsorted(bounds[i], ct_lower_bounds[i]) - 1 for i in range(3)])
    lower_indices[lower_indices < 0] = 0

    upper_indices = np.array([np.searchsorted(bounds[i], ct_upper_bounds[i]) - 1 for i in range(3)])
    upper_indices[upper_indices < 0] = 0

    return lower_indices, upper_indices


def calculate_weights_of_ct_voxel(voxel_ijk)