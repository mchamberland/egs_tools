import operator
import numpy as np
from typing import List, Tuple


def get_voxel_center_from_index(index: int, bounds3d: List[List[float]]):
    """Return the coordinates of the center of the voxel specified with its 1-D index

    Args:
        index (int): 1-D index
        bounds3d (list[list[float]]): list of list of the bounds of each voxel in each dimension

    Returns:
        (float, float, float): coordinates of the center of the voxel
    """
    dimensions = [len(bounds3d[0]) - 1, len(bounds3d[1]) - 1, len(bounds3d[2]) - 1]
    return tuple((bounds3d[d][i] + bounds3d[d][i + 1]) / 2.
                 for d, i in enumerate(get_ijk_from_index(index, dimensions)))


def get_voxel_center_from_ijk(index3d: Tuple[int, int, int], bounds3d: List[List[float]]):
    """Return the coordinates of the center of the voxel specified with its multi-dimensional index

    Args:
        index3d (Tuple(int)): multi-dimensional index
        bounds3d (list[list[float]]): list of list of the bounds of each voxel in each dimension

    Returns:
        (Tuple(float)): coordinates (x,y,z) of the center of the voxel
    """
    return tuple((bounds3d[d][i] + bounds3d[d][i + 1]) / 2. for d, i in enumerate(index3d))


def get_pixel_center_from_ij(index2d: Tuple[int, int], bounds2d: List[List[float]]):
    """Return the coordinates of the center of the pixel specified with its multi-dimensional index

    Args:
        index2d (Tuple(int)): multi-dimensional index
        bounds2d (list[list[float]]): list of list of the bounds of each pixel in each dimension

    Returns:
        (Tuple(float)): coordinates (x,y) of the center of the pixel
    """
    return tuple((bounds2d[d][i] + bounds2d[d][i + 1]) / 2. for d, i in enumerate(index2d))


def get_all_pixel_centers(bounds2d: List[List[float]]):
    """Return a list of the coordinates of the center of all pixels

    Args:
        bounds2d (list[list[float]]): list of list of the bounds of each pixel in each dimension

    Returns:
        list[tuple[float]]: list of coordinates (x,y) of the center of every pixel
    """
    return [get_pixel_center_from_ij((i, j), bounds2d) for j in range(len(bounds2d[1]) - 1)
            for i in range(len(bounds2d[0]) - 1)]


def get_voxel_size_from_ijk(index3d: Tuple[int, int, int], bounds3d: List[List[float]]):
    """Return the size of the rectilinear voxel specified by its multi-dimensional index

    Args:
        index3d (Tuple(int)): multi-dimensional index
        bounds3d (list[list[float]]): list of list of the bounds of each voxel in each dimension

    Returns:
        (Tuple(float)): size of the voxel, in each dimension
    """
    return tuple(bounds3d[d][i + 1] - bounds3d[d][i] for d, i in enumerate(index3d))


def get_index_from_ijk(index3d: Tuple[int, int, int], dimensions: List[int]):
    """Return the 1-D index corresponding to the (i,j,k) index

    Args:
        index3d (Tuple(int)): multi-dimensional index
        dimensions (list[int]): the dimensions of the voxel grid
    Returns:
        (int): 1-D index of voxel
    """
    i, j, k = index3d
    return i + j * dimensions[0] + k * dimensions[0] * dimensions[1]


def get_ijk_from_xyz(pos: Tuple[float, float, float], bounds3d: List[List[float]]):
    """Return the (i, j, k) index of the voxel containing position (x, y, z)

    Args:
        pos (tuple(float)): (x, y, z) position
        bounds3d (list[list[float]]): list of list of the bounds of each voxel in each dimension

    Returns:
        (tuple(int)): multi-dimensional index (i, j, k) of voxel
    """
    return tuple(get_index_from_position(p, bounds3d[i]) for i, p in enumerate(pos))


def get_ij_from_xy(pos: Tuple[float, float], bounds2d: List[List[float]]):
    """Return the (i, j) index of the pixel containing position (x, y)

    Args:
        pos (tuple(float)): (x, y) position
        bounds2d (list[list[float]]): list of list of the bounds of each pixel in each dimension

    Returns:
        (tuple(int)): multi-dimensional index (i, j) of voxel
    """
    return tuple(get_index_from_position(p, bounds2d[i]) for i, p in enumerate(pos))


def get_ijk_from_index(index: int, dimensions: List[int]):
    """Return the i,j,k index of the 1-D voxel index

    Args:
        index (int): 1-D voxel index
        dimensions (list[int]): the dimensions of the voxel grid

    Returns:
        (int, int, int): i,j,k index of the voxel
    """

    i = operator.mod(index, dimensions[0])
    ij = dimensions[0] * dimensions[1]
    k = (index - i) // ij
    j = (index - i - k * ij) // dimensions[0]

    return i, j, k


def get_ijk_indexing_array_from_index_list(index_list: List[int], dimensions: List[int]):
    """Return an ndarray that can be used to index a multidimensional array to access the elements pointed to by the
    list of 1-D indices

    Args:
        index_list (numpy.ndarray[int]): list of 1-D voxel index
        dimensions (list[int]): the dimensions of the voxel grid

    Returns:
        (list[list[int]]): list of 3 lists of integers that allows to index a multidimensional array
    """

    voxels = [get_ijk_from_index(v, dimensions) for v in index_list]

    return list(map(list, zip(*voxels)))


def get_index_from_position(pos: float, bounds: List[float], is_strictly_inside=True):
    """Return the 1-D index of the voxel containing the position

    If is_strictly_inside is set to true (default), then a position outside the bounds will return -1.

    If is_strictly_inside is set to false, then a position below the lower bound will return the lowest index and
    a position above the upper bound will return the highest bound.

    Args:
        pos (float): 1-D position
        bounds (list[float]): bounds of the voxels
        is_strictly_inside (bool): index returned only if position is strictly inside

    Returns:
        (int): 1-D index of the voxel which contains the position
    """
    # outside all bounds
    if is_strictly_inside and ((pos < bounds[0]) or (bounds[-1] < pos)):
        return -1
    elif not is_strictly_inside:
        if pos < bounds[0]:
            return 0
        elif bounds[-1] < pos:
            return len(bounds) - 2

    # check end cases
    if pos == bounds[0]:
        return 0
    elif pos == bounds[-1]:
        return len(bounds) - 2

    # normal case
    idx = 0
    while pos > bounds[idx]:
        idx += 1

    return idx - 1


def get_region_indices_from_low_and_hi_pos(low_pos: float, high_pos: float, bounds: List[float],
                                           include_partial_voxel=True):
    """Return a tuple of the lower and higher 1-D indices
    of the voxels defining the region of the dose distribution.

    In other words, this tells you which voxels define the bounds (low_pos to high_pos)
    of a region of the dose distribution.

    If include_partial_voxel is set to false, then a position inside a voxel will return the nearest innermost
    bound.

    If include_partial_voxel is set to true (default), then a position inside a voxel will return the index of this
    voxel.

    Args:
        low_pos (float): 1-D position of the lower bound of the region
        high_pos (float): 1-D position of the higher bound of the region
        bounds (list[float]): bounds of the voxels
        include_partial_voxel (bool): include partial voxels in the region or not

    Returns:
        (int, int): 1-D indices of the lower and higher bounds of the region
    """
    # outside all bounds or wrong low/high pos
    if bounds[-1] < low_pos or high_pos < bounds[0] or low_pos > high_pos:
        return -1, -1

    # check high end cases
    if high_pos > bounds[-1]:
        higher_bound_idx = get_index_from_position(bounds[-1], bounds, is_strictly_inside=False)
    else:
        higher_bound_idx = get_index_from_position(high_pos, bounds)
        if high_pos == bounds[higher_bound_idx] or include_partial_voxel:
            pass
        elif not include_partial_voxel and high_pos < bounds[-1]:
            higher_bound_idx -= 1

    # check lower end cases
    if low_pos <= bounds[0]:
        return 0, higher_bound_idx
    elif low_pos == bounds[-1]:
        return len(bounds) - 2, higher_bound_idx

    lower_bound_idx = 0
    while low_pos >= bounds[lower_bound_idx]:
        lower_bound_idx += 1

    if include_partial_voxel:
        return lower_bound_idx - 1, higher_bound_idx
    else:
        return lower_bound_idx, higher_bound_idx


def get_next_x_neighbour(index3d: Tuple[int, int, int], x_bounds: List[float]):
    """Return the 1-D index of the next x neighbour of the voxel specified by the (i,j,k) index

    Args:
        index3d (Tuple(int)): multi-dimensional (i,j,k) index of voxel
        x_bounds (list[float]): list of the bounds of the voxels in x

    Returns:
        (int): 1-D index of next x neighbour
    """
    i, j, k = index3d

    if i > len(x_bounds) - 2:
        return -1
    elif i == len(x_bounds) - 2:
        return i
    else:
        return i + 1


def get_prev_x_neighbour(index3d: Tuple[int, int, int]):
    """Return the 1-D index of the previous x neighbour of the voxel specified by the index

    Args:
        index3d (Tuple(int)): multi-dimensional (i,j,k) index of voxel

    Returns:
        (int): 1-D index of previous x neighbour
    """
    i, j, k = index3d

    if i == 0:
        return i
    else:
        return i - 1


def get_next_y_neighbour(index3d: Tuple[int, int, int], y_bounds: List[float]):
    """Return the 1-D index of the next y neighbour of the voxel specified by the index

    Args:
        index3d (Tuple(int)): multi-dimensional (i,j,k) index of voxel
        y_bounds (list[float]): list of the bounds of the voxels in y

    Returns:
        (int): 1-D index of next y neighbour
    """
    i, j, k = index3d

    if j > len(y_bounds) - 2:
        return -1
    elif j == len(y_bounds) - 2:
        return j
    else:
        return j + 1


def get_prev_y_neighbour(index3d: Tuple[int, int, int]):
    """Return the 1-D index of the previous y neighbour of the voxel specified by the index

    Args:
        index3d (Tuple(int)): multi-dimensional (i,j,k) index of voxel

    Returns:
        (int): 1-D index of previous y neighbour
    """
    i, j, k = index3d

    if j == 0:
        return j
    else:
        return j - 1


def get_next_z_neighbour(index3d: Tuple[int, int, int], z_bounds: List[float]):
    """Return the 1-D index of the next z neighbour of the voxel specified by the index

    Args:
        index3d (Tuple(int)): multi-dimensional (i,j,k) index of voxel
        z_bounds (list[float]): list of the bounds of the voxels in z

    Returns:
        (int): 1-D index of next z neighbour
    """
    i, j, k = index3d

    if k > len(z_bounds) - 2:
        return -1
    elif k == len(z_bounds) - 2:
        return k
    else:
        return k + 1


def get_prev_z_neighbour(index3d: Tuple[int, int, int]):
    """Return the 1-D index of the previous z neighbour of the voxel specified by the index

    Args:
        index3d (Tuple(int)): multi-dimensional (i,j,k) index of voxel

    Returns:
        (int): 1-D index of previous z neighbour
    """
    i, j, k = index3d

    if k == 0:
        return k
    else:
        return k - 1


def are_bounds_within_tolerance(bounds3d1: List[List[float]], bounds3d2: List[List[float]]) -> bool:
    return np.all([np.allclose(bounds3d1[i], bounds3d2[i]) for i in range(len(bounds3d1))])
