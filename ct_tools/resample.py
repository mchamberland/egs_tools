import numpy as np
from typing import Tuple
from ct_tools.ctdata import CTdata


def resample_ctdata(ct: CTdata, voxels: Tuple[float, float, float], size_or_count='size') -> CTdata:
    # based on the the resampleCT subroutine in ctcreate.mortran in EGSnrc
    # resampling weighs the original CT values by the fractional volume of each resampled (xyz) voxel that overlaps
    # with a given CT voxel
    if size_or_count == "size":
        voxel_size_in_cm = voxels
    elif size_or_count == "count":
        voxel_size_in_cm = [ct.image_size_in_cm()[i] / v for (i, v) in enumerate(voxels)]
    else:
        raise Exception("Argument 'size_or_count' must be either 'size' or 'count'.")
    print("Resampling CT data...")
    print("Requested voxel size in cm: ({:.4f}, {:.4f}, {:.4f})".format(voxel_size_in_cm[0],
                                                                        voxel_size_in_cm[1],
                                                                        voxel_size_in_cm[2]))

    check_requested_voxel_size(ct.image_size_in_cm(), voxel_size_in_cm)
    dimensions, adjusted_voxel_size = adjust_requested_voxel_size(ct.image_size_in_cm(), voxel_size_in_cm)
    bounds = calculate_new_bounds(ct.bounds, adjusted_voxel_size, dimensions)

    resampled = CTdata()
    resampled.dimensions = dimensions
    resampled.bounds = bounds
    resampled.image = np.zeros(dimensions)

    print("Calculating weights and new CT values. This may take a while...")
    loop_counter = 0
    print_counter = 1
    for (ct_ijk, value) in np.ndenumerate(ct.image):
        loop_counter += 1
        if loop_counter % int(ct.nvox() / 10) == 0:
            print("{:d}%...".format(print_counter * 10))
            print_counter += 1
        lower_ijk, upper_ijk = find_ijk_where_ct_boundaries_lie(ct.bounds, bounds, ct_ijk)
        indices = list(zip(lower_ijk, upper_ijk))
        weights = calculate_weights_of_ct_voxel(ct.bounds, bounds, ct_ijk, indices, ct.voxel_size_in_cm(),
                                                adjusted_voxel_size)

        for i in range(lower_ijk[0], upper_ijk[0] + 1):
            for j in range(lower_ijk[1], upper_ijk[1] + 1):
                for k in range(lower_ijk[2], upper_ijk[2] + 1):
                    resampled.image[(i, j, k)] += ct.image[ct_ijk] * weights[0][i] * weights[1][j] * weights[2][k]

    print("Resampling completed!")
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
    print("New voxel size in cm: ({:.4f}, {:.4f}, {:.4f})".format(adjusted_voxel_size[0],
                                                                  adjusted_voxel_size[1],
                                                                  adjusted_voxel_size[2]))
    print("New number of voxels: ({:d}, {:d}, {:d})".format(dimensions[0],
                                                            dimensions[1],
                                                            dimensions[2]))
    return dimensions, adjusted_voxel_size


def calculate_new_bounds(ct_bounds, voxel_size, dimensions):
    print("Calculating new bounds...")
    return [[ct_bounds[i][0] + j * voxel_size[i] for j in range(dimensions[i] + 1)] for i in range(3)]


def find_ijk_where_ct_boundaries_lie(ct_bounds, bounds, ct_ijk):
    ct_lower_bounds = [ct_bounds[i][v] for (i, v) in enumerate(ct_ijk)]
    ct_upper_bounds = [ct_bounds[i][v + 1] for (i, v) in enumerate(ct_ijk)]

    lower_ijk = np.array([np.searchsorted(bounds[i], ct_lower_bounds[i]) - 1 for i in range(3)])
    lower_ijk[lower_ijk < 0] = 0

    upper_ijk = np.array([np.searchsorted(bounds[i], ct_upper_bounds[i]) - 1 for i in range(3)])
    upper_ijk[upper_ijk < 0] = 0

    return lower_ijk, upper_ijk


def calculate_weights_of_ct_voxel(ct_bounds, bounds, ct_ijk, indices, ct_voxel_size, voxel_size):
    max_xindex = max(len(ct_bounds[0]), len(bounds[0]))
    max_yindex = max(len(ct_bounds[1]), len(bounds[1]))
    max_zindex = max(len(ct_bounds[2]), len(bounds[2]))

    weights = np.array([np.zeros(max_xindex), np.zeros(max_yindex), np.zeros(max_zindex)])

    for (dimension, index) in enumerate(indices):
        if index[0] == index[1]:
            # ct low and hi bounds are in same xyz voxel
            weights[dimension][index[0]] = ct_voxel_size[dimension] / voxel_size[dimension]
        else:
            for i in range(index[0], index[1] + 1):
                if bounds[dimension][i] >= ct_bounds[dimension][ct_ijk[dimension]] and \
                                bounds[dimension][i + 1] <= ct_bounds[dimension][ct_ijk[dimension] + 1]:
                    # xyz voxel is entirely in ct voxel
                    weights[dimension][i] = 1.0

                elif bounds[dimension][i] <= ct_bounds[dimension][ct_ijk[dimension]] and \
                                bounds[dimension][i + 1] <= ct_bounds[dimension][ct_ijk[dimension] + 1]:
                    # ct voxel straddles the upper bound of the xyz voxel
                    weights[dimension][i] = (bounds[dimension][i + 1] - ct_bounds[dimension][ct_ijk[dimension]]) / \
                                            voxel_size[dimension]

                else:
                    # ct voxel straddles the lower bound of the xyz voxel
                    weights[dimension][i] = (ct_bounds[dimension][ct_ijk[dimension] + 1] - bounds[dimension][i]) / \
                                            voxel_size[dimension]

    return weights
