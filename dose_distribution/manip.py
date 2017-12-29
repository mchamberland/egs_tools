import math
import voxelnav
import numpy as np
from typing import Tuple
from .dose3d import DoseDistribution


def average_about_x(the_dose3d: DoseDistribution):
    """Average the dose in voxels symmetric about the x = 0 plane"""
    nx, ny, nz = the_dose3d.dimensions
    num_steps = nx // 2

    for jj in range(ny):
        for kk in range(nz):
            for ii in range(num_steps):
                min_index = (ii, jj, kk)
                max_index = (nx - ii - 1, jj, kk)
                dose, unc = average_dose(the_dose3d, min_index, max_index)
                the_dose3d.dose[min_index], the_dose3d.fract_unc[min_index] = dose, unc
                the_dose3d.dose[max_index], the_dose3d.fract_unc[max_index] = dose, unc


def average_about_y(the_dose3d: DoseDistribution):
    """Average the dose in voxels symmetric about the y = 0 plane"""
    nx, ny, nz = the_dose3d.dimensions
    num_steps = ny // 2

    for ii in range(nx):
        for kk in range(nz):
            for jj in range(num_steps):
                min_index = (ii, jj, kk)
                max_index = (ii, ny - jj - 1, kk)
                dose, unc = average_dose(the_dose3d, min_index, max_index)
                the_dose3d.dose[min_index], the_dose3d.fract_unc[min_index] = dose, unc
                the_dose3d.dose[max_index], the_dose3d.fract_unc[max_index] = dose, unc


def average_about_z(the_dose3d: DoseDistribution):
    """Average the dose in voxels symmetric about the z = 0 plane in a rectilinear dose distribution"""
    nx, ny, nz = the_dose3d.dimensions
    num_steps = nz // 2

    for ii in range(nx):
        for jj in range(ny):
            for kk in range(num_steps):
                min_index = (ii, jj, kk)
                max_index = (ii, jj, nz - kk - 1)
                dose, unc = average_dose(the_dose3d, min_index, max_index)
                the_dose3d.dose[min_index], the_dose3d.fract_unc[min_index] = dose, unc
                the_dose3d.dose[max_index], the_dose3d.fract_unc[max_index] = dose, unc


def average_dose(the_dose3d: DoseDistribution, idx_1: Tuple[int, int, int], idx_2: Tuple[int, int, int]):
    """Return the average of the dose and its fractional uncertainty in voxels specified by two indexes

    This returns the weighted average of the two doses, with weights equal to (1 / fractional uncertainty)^2,
    and fractional uncertainty equal to (1 / sqrt(sum of the weights)).

    Two special cases:
        Return 0, 0 if any of the two doses are equal to 0.
        Return (dose_1 + dose_2) / 2 with 0 uncertainty if any of the uncertainties are equal to 0.

    Args:
        the_dose3d (DoseDistribution): a 3ddose instance
        idx_1 (Tuple(int)): multi-dimensional index of first voxel
        idx_2 (Tuple(int)): multi_dimensional index of second voxel

    Returns:
        (float, float): the average dose and its fractional uncertainty
    """
    dose_values = np.array([the_dose3d.dose[idx_1], the_dose3d.dose[idx_2]])
    unc_values = np.array([the_dose3d.fract_unc[idx_1], the_dose3d.fract_unc[idx_2]])

    if not dose_values.all():
        return 0, 0

    if not unc_values.all():
        return dose_values.mean(), 0
    else:
        weights = (1 / unc_values) ** 2

    return calculate_weighted_average(dose_values, weights)


def calculate_weighted_average(values: np.ndarray, weights: np.ndarray):
    """Return the weighted average of the values and its uncertainty

    Calculate the weighted average of the values held in 'values', with weights according
    to 'weights'. The uncertainty is given as (1 / sqrt(sum of the weights)).

    Args:
        values (numpy.ndarray): values to be averaged
        weights (numpy.ndarray): weights for each value

    Returns
        (float, float): the weighted average and its uncertainty
    """
    the_average, sum_of_weights = np.average(values, weights=weights, returned=True)
    the_uncertainty = 1 / math.sqrt(sum_of_weights)

    return the_average, the_uncertainty


def tri_linear_interp_dose_xyz(the_dose3d: DoseDistribution, pos: Tuple[float, float, float],
                               interp_fun=lambda x, y, z: 1):
    """Return the tri-linearly interpolated dose at position (x, y, z)

    Uses the 8 nearest neighbouring voxels to interpolate the dose at (x, y, z). Before interpolation,
    all doses are multiplied by the interpolation function interp_fun. This method is only used for
    rectilinear dose distributions.

    If the interpolated dose is equal to 0, the fractional uncertainty is set to 1.

    The uncertainty on the interpolated dose is taken as the uncertainty on a weighted arithmetic mean.
    See here for more details: https://en.wikipedia.org/wiki/Weighted_arithmetic_mean

    Args:
        the_dose3d (DoseDistribution): a 3ddose instance
        pos (tuple(float)): the point (x,y,z) of interest
        interp_fun (function = x,y,z = 1): before interpolation, the neighbouring doses
        will be multiplied by this function

    Returns:
        (float, float): the interpolated dose and its fractional uncertainty
    """
    # 1-D indices of nearest neighbours
    idx = voxelnav.get_ijk_from_xyz(pos, the_dose3d.bounds)
    ii, jj, kk = idx
    x, y, z = pos
    vox_x, vox_y, vox_z = voxelnav.get_voxel_center_from_ijk(idx, the_dose3d.bounds)
    if x < vox_x:
        i_min = voxelnav.get_prev_x_neighbour(idx)
        i_max = ii
    else:
        i_min = ii
        i_max = voxelnav.get_next_x_neighbour(idx, the_dose3d.bounds[0])

    if y < vox_y:
        j_min = voxelnav.get_prev_y_neighbour(idx)
        j_max = jj
    else:
        j_min = jj
        j_max = voxelnav.get_next_y_neighbour(idx, the_dose3d.bounds[1])

    if z < vox_z:
        k_min = voxelnav.get_prev_z_neighbour(idx)
        k_max = kk
    else:
        k_min = kk
        k_max = voxelnav.get_next_z_neighbour(idx, the_dose3d.bounds[2])

    index_list = [(i_min, j_min, k_min),
                  (i_max, j_min, k_min),
                  (i_min, j_max, k_min),
                  (i_max, j_max, k_min),
                  (i_min, j_min, k_max),
                  (i_max, j_min, k_max),
                  (i_min, j_max, k_max),
                  (i_max, j_max, k_max)]

    dose_values = np.zeros(len(index_list))
    abs_unc_squared = np.zeros(len(index_list))
    for i, index in enumerate(index_list):
        xint, yint, zint = voxelnav.get_voxel_center_from_ijk(index, the_dose3d.bounds)
        dose_values[i] = interp_fun(xint, yint, zint) * the_dose3d.dose[index]
        abs_unc = the_dose3d.fract_unc[index] * the_dose3d.dose[index]
        abs_unc_squared[i] = abs_unc ** 2

    # coordinates of centres of nearest neighbours
    x_min, y_min, z_min = voxelnav.get_voxel_center_from_ijk((i_min, j_min, k_min), the_dose3d.bounds)
    x_max, y_max, z_max = voxelnav.get_voxel_center_from_ijk((i_max, j_max, k_max), the_dose3d.bounds)

    max_x_diff = x_max - x
    min_x_diff = x - x_min
    max_y_diff = y_max - y
    min_y_diff = y - y_min
    max_z_diff = z_max - z
    min_z_diff = z - z_min

    weights = np.array([max_x_diff * max_y_diff * max_z_diff,
                        min_x_diff * max_y_diff * max_z_diff,
                        max_x_diff * min_y_diff * max_z_diff,
                        min_x_diff * min_y_diff * max_z_diff,
                        max_x_diff * max_y_diff * min_z_diff,
                        min_x_diff * max_y_diff * min_z_diff,
                        max_x_diff * min_y_diff * min_z_diff,
                        min_x_diff * min_y_diff * min_z_diff])

    weights_squared = weights ** 2

    dose = np.average(dose_values, weights=weights) / interp_fun(x, y, z)
    abs_unc = math.sqrt(np.average(abs_unc_squared, weights=weights_squared))

    if dose != 0:
        return dose, abs_unc / dose
    else:
        return 0, 1

# TODO add dose comparison methods, similar to 3ddose_tools.
# TODO add histogram plotting for comparison methods
# TODO add pygrace stuff to produce xmgrace plots
