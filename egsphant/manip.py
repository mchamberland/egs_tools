import voxelnav
from .egsphant import EGSphant


def trim_phantom(the_egsphant: EGSphant, the_bounds):
    """Trim the egsphant to the bounds specified

    Args:
        the_egsphant (EGSphant): the ESGphant instance to trim
        the_bounds (tuple(float)): the bounds of the region of interest

    """
    if len(the_bounds) > 1:
        xi, xf, yi, yf, zi, zf = the_bounds
    else:
        xi = yi = zi = (-1) * the_bounds[0] / 2.
        xf = yf = zf = the_bounds[0] / 2.

    x_index = voxelnav.get_region_indices_from_low_and_hi_pos(xi, xf, the_egsphant.bounds[0])
    y_index = voxelnav.get_region_indices_from_low_and_hi_pos(yi, yf, the_egsphant.bounds[1])
    z_index = voxelnav.get_region_indices_from_low_and_hi_pos(zi, zf, the_egsphant.bounds[2])

    xb_slice = slice(x_index[0], x_index[1] + 2)
    yb_slice = slice(y_index[0], y_index[1] + 2)
    zb_slice = slice(z_index[0], z_index[1] + 2)

    the_egsphant.bounds[0] = the_egsphant.bounds[0][xb_slice]
    the_egsphant.bounds[1] = the_egsphant.bounds[1][yb_slice]
    the_egsphant.bounds[2] = the_egsphant.bounds[2][zb_slice]

    x_slice = slice(x_index[0], x_index[1] + 1)
    y_slice = slice(y_index[0], y_index[1] + 1)
    z_slice = slice(z_index[0], z_index[1] + 1)

    the_egsphant.phantom = the_egsphant.phantom[x_slice, y_slice, z_slice]
    the_egsphant.density = the_egsphant.density[x_slice, y_slice, z_slice]

    the_egsphant.dimensions[0] = len(the_egsphant.bounds[0]) - 1
    the_egsphant.dimensions[1] = len(the_egsphant.bounds[1]) - 1
    the_egsphant.dimensions[2] = len(the_egsphant.bounds[2]) - 1


def copy_medium_from_source_to_destination(medium: str, source: EGSphant, destination: EGSphant):
    """All instances of 'medium' in 'source' will replace the medium in the corresponding voxels in 'destination'

    If the two egsphants have different dimensions/voxel sizes, then voxels are mapped by the location of their centre.

    Args:
        medium (str): the medium that will replace the current instance
        source (EGSphant): the source egsphant containing the medium to be copied
        destination (EGSphant): the destination egsphant to modify

    """
    voxels_with_medium = source.get_index_of_voxels_with_medium(medium)
    voxels = voxelnav.get_ijk_indexing_array_from_index_list(voxels_with_medium, source.dimensions)

    density = source.density[voxels]
    medium_key = add_medium(destination, medium)

    if source.bounds == destination.bounds:
        destination.phantom[voxels] = medium_key
        destination.density[voxels] = density
    else:
        for v in voxels_with_medium:
            pos = voxelnav.get_voxel_center_from_index(v, source.bounds)
            ijk = voxelnav.get_ijk_from_xyz(pos, destination.bounds)
            density = source.density[voxelnav.get_ijk_from_index(v, source.dimensions)]
            destination.phantom[ijk] = medium_key
            destination.density[ijk] = density


def replace_original_medium_with_new_medium(the_egsphant: EGSphant, original_medium: str, new_medium: str,
                                            density: float):
    """All instances of 'original_medium' will be replaced by 'medium' with nominal density 'density'

    Args:
        the_egsphant (EGSphant): the egsphant to modify
        original_medium (str): the original medium to replace
        new_medium (str): the medium that will replace the original medium
        density (float): nominal density of medium in g / cm^3

    """
    # Get the key of the original medium, the indices of the voxels to replace, and the density of the medium
    original_medium_key = the_egsphant.get_medium_key(original_medium)
    voxels_to_replace = [i for i, x in enumerate(the_egsphant.phantom.flatten(order='F')) if x == original_medium_key]

    # Figure out if the medium already exists in this egsphant or if it needs to be added
    index = [i for i, m in enumerate(the_egsphant.media) if m == new_medium]
    if len(index) > 0:
        medium_key = the_egsphant.get_medium_key(new_medium)
    else:
        medium_key = the_egsphant.MEDIUM_KEY_STRING[the_egsphant.number_of_media]
        the_egsphant.number_of_media += 1
        the_egsphant.media.append(new_medium)
        the_egsphant.medium_keys[medium_key] = new_medium
        the_egsphant.inverse_key_mapping[new_medium] = medium_key

    voxels = voxelnav.get_ijk_indexing_array_from_index_list(voxels_to_replace, the_egsphant.dimensions)
    the_egsphant.phantom[voxels] = medium_key
    the_egsphant.density[voxels] = density


def add_medium(the_egsphant: EGSphant, medium: str):
    """Add the medium to the list of media of the egsphant and update the number of media and the medium key mapping.

    If the medium already exists in the egsphant, nothing is done. In both cases, the medium key is returned.

    Args:
        the_egsphant (EGSphant): the egsphant instance
        medium (str): the medium that will be added (if needed)

    Returns:
        str: medium key

    """
    index = [i for i, m in enumerate(the_egsphant.media) if m == medium]
    if len(index) > 0:
        medium_key = the_egsphant.get_medium_key(medium)
    else:
        medium_key = the_egsphant.MEDIUM_KEY_STRING[the_egsphant.number_of_media]
        the_egsphant.number_of_media += 1
        the_egsphant.media.append(medium)
        the_egsphant.medium_keys[medium_key] = medium
        the_egsphant.inverse_key_mapping[medium] = medium_key

    return medium_key
