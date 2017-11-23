from typing import Tuple
from ct_tools.ctdata import CTdata


def resample_ctdata(original: CTdata, new_voxel_size_in_cm: Tuple[float, float, float]) -> CTdata:
    resampled = CTdata()

    image_size = original.image_size_in_cm()
    if any(v < 0 for v in new_voxel_size_in_cm):
        raise Exception("Voxel size must be positive.")

    if any(v > image_size[i] for (i, v) in enumerate(new_voxel_size_in_cm)):
        raise Exception("Voxel size in any one direction cannot be greater"
                        "than the size of the original CT in this direction.")

    dimensions_new = [int(image_size[i] / v) for (i, v) in enumerate(new_voxel_size_in_cm)]
    adjusted_voxel_size = [float(image_size[i] / n) for (i, n) in enumerate(dimensions_new)]
    print("Adjusted dimensions (in cm) so that an integer number of voxels fit exactly on the CT data:")
    print("New voxel size in cm: ({:.2f}, {:.2f}, {:.2f})".format(adjusted_voxel_size[0],
                                                                  adjusted_voxel_size[1],
                                                                  adjusted_voxel_size[2]))
    print("New number of voxels: ({:d}, {:d}, {:d})".format(dimensions_new[0],
                                                            dimensions_new[1],
                                                            dimensions_new[2]))

    bounds_new = [[original.bounds[i][0] + j * adjusted_voxel_size[i] for j in range(dimensions_new[i] + 1)]
                  for i in range(3)]

    dimensions = original.dimensions
    voxel_size = original.voxel_size_in_cm()
    bounds = original.bounds

    return resampled
