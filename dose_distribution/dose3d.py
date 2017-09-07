"""Dose distribution

This module exports a single class that handles dose distributions from 3ddose files.

class:

DoseDistribution -- dose distribution from a 3ddose file
"""
import os
import gzip
from typing import Tuple
import numpy as np
import voxelnav


class DoseDistribution:
    """Dose distribution from a 3ddose file

    The class stores a dose distribution from a 3ddose file whose filename is passed as an argument.
    The expected format of the 3ddose file is as follows:

        1) 3 integers representing the number of voxels in each dimension. For a rectilinear distribution,
        this corresponds to nx, ny, nz.

        2) The bounds of the voxels in each dimension.

        3) The dose in each voxel in Gy per hist, using Fortran-style indexing where the index of the first dimension
        changes the fastest.

        4) The fractional uncertainty on the dose in each voxel.

    Attributes:
        dimensions (list[int]): number of voxels in each dimension of the distribution
        bounds (list[list[float]]): list of list of the bounds of each voxel in each dimension
        dose (numpy.ndarray): multi-dimensional (nx x ny x nz; nz x nr; or nr) array storing the dose
        fract_unc (numpy.ndarray): multi-dimensional array storing the fractional uncertainty on the dose

    Methods:
        __init__(filename)
            Create a DoseDistribution instance from a 3ddose filename

        read_dose(filename)
            Read a 3ddose file specified by its filename

        write_dose(filename)
            Write a dose distribution to a 3ddose file specified by filename

        nvox()
            Return the total number of voxels

        get_max_dose()
            Return the maximum dose of the distribution and its fractional uncertainty

        get_min_dose()
            Return the minimum dose of the distribution and its fractional uncertainty

        get_point_dose(pos)
            Return the dose and its fractional uncertainty in the voxel containing position (x, y, z)

        get_region_dose(the_bounds)
            Return the dose and its fractional uncertainty in the region of the dose distribution defined
            by the bounds

    """

    def __init__(self, filename=None):
        """Create a DoseDistribution instance from a 3ddose filename"""
        self.dimensions = []
        self.bounds = []
        self.dose = np.array([0])
        self.fract_unc = np.array([0])

        if filename is not None:
            self.dose, self.fract_unc = self.read_dose(filename)

    def read_dose(self, filename) -> (np.ndarray, np.ndarray):
        """Read a txt or gzip 3ddose file specified by its filename

        This method fills out all the attributes of the DoseDistribution instance.
        After calling this method, the dose and fractional uncertainty arrays of the DoseDistribution instance
        are stored in numpy arrays with dimensions (nx x ny x nz).

        Args:
            filename (str): filename of the 3ddose file to be read

        Returns
            (numpy.ndarray, numpy.ndarray): dose and its fractional uncertainty
        """
        if not os.path.exists(filename):
            return -1

        if filename.endswith('gz'):
            with gzip.open(filename, 'r') as file:
                lines = file.readlines()
        else:
            with open(filename, 'r') as file:
                lines = file.readlines()

        self.dimensions = list(map(int, lines[0].split()))
        self.bounds = [list(map(float, line.split())) for line in lines[1:4]]

        dose = np.array(list(map(float, lines[4].split())))
        dose = dose.reshape(self.dimensions, order='F')

        fract_unc = np.array(list(map(float, lines[5].split())))
        fract_unc = fract_unc.reshape(self.dimensions, order='F')

        return dose, fract_unc

    def write_dose(self, filename):
        """Write a dose distribution to a 3ddose file specified by filename

        Args:
            filename (str): filename of the 3ddose file to be written
        """

        if filename.endswith('gz'):
            file = gzip.open(filename, 'wt')
        else:
            file = open(filename, 'wt')

        file.write(' '.join(str(vox) for vox in self.dimensions) + '\n')
        file.write('\n'.join(' '.join(str(b) for b in d) for d in self.bounds) + '\n')
        file.write(' '.join('{:.8f}'.format(d) for d in self.dose.flatten('F')) + '\n')
        file.write(' '.join('{:.8f}'.format(u) for u in self.fract_unc.flatten('F')) + '\n')
        file.close()

    def nvox(self) -> int:
        """Return the total number of voxels"""
        nvox = 1
        for dims in self.dimensions:
            nvox *= dims
        return nvox

    def get_max_dose(self) -> (float, float):
        """Return the maximum dose of the distribution and its fractional uncertainty

        Returns
            (float, float): the maximum dose and its fractional uncertainty
        """
        max_index = np.unravel_index(np.argmax(self.dose), tuple(self.dimensions))

        return np.max(self.dose), self.fract_unc[max_index]

    def get_min_dose(self) -> (float, float):
        """Return the minimum dose of the distribution and its fractional uncertainty

        Returns
            (float, float): the minimum dose and its fractional uncertainty
        """
        min_index = np.unravel_index(np.argmin(self.dose), tuple(self.dimensions))

        return np.min(self.dose), self.fract_unc[min_index]

    def get_point_dose(self, pos: Tuple[float, float, float]) -> (float, float):
        """Return the dose and its fractional uncertainty in the voxel containing position (x, y, z)

        Args:
            pos (tuple(float)): the position (x, y, z)

        Returns:
            (float, float): the dose and its fractional uncertainty
        """
        index = voxelnav.get_ijk_from_xyz(pos, self.bounds)

        return self.dose[index], self.fract_unc[index]

    def get_region_dose(self, the_bounds) -> (np.ndarray, np.ndarray):
        """Return the region of the dose and its fractional uncertainty in
        the region of the dose distribution defined by the bounds

        Args:
            the_bounds (tuple(float)): the bounds of the region of interest

        Returns:
            (numpy.ndarray, numpy.ndarray): the dose and its fractional uncertainty in the region
        """
        if len(the_bounds) > 1:
            xi, xf, yi, yf, zi, zf = the_bounds
        else:
            xi = yi = zi = (-1) * the_bounds[0] / 2.
            xf = yf = zf = the_bounds[0] / 2.

        x_index = voxelnav.get_region_indices_from_low_and_hi_pos(xi, xf, self.bounds[0])
        y_index = voxelnav.get_region_indices_from_low_and_hi_pos(yi, yf, self.bounds[1])
        z_index = voxelnav.get_region_indices_from_low_and_hi_pos(zi, zf, self.bounds[2])

        x_slice = slice(x_index[0], x_index[1] + 1)
        y_slice = slice(y_index[0], y_index[1] + 1)
        z_slice = slice(z_index[0], z_index[1] + 1)

        return self.dose[x_slice, y_slice, z_slice], self.fract_unc[x_slice, y_slice, z_slice]
