"""Dose distribution

This module exports a single class that handles dose distributions from 3ddose files.

class:

DoseDistribution -- dose distribution from a 3ddose file
"""
import os
import gzip
import pydicom
from os.path import join
from typing import Tuple
import numpy as np
import voxelnav
import value_mapping as vmap
import dicom.uid
import random

EGS_TOOLS_HOME = os.path.expandvars("$EGS_TOOLS_HOME")
EMPTY_DICOM_TEMPLATE = "RT_Dose_template.dcm"


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
                values = [float(v) for line in file.readlines() for v in line.split()]
        else:
            with open(filename, 'r') as file:
                values = [float(v) for line in file.readlines() for v in line.split()]

        self.dimensions = [int(values[i]) for i in range(0, 3)]

        xmin = 3
        xmax = xmin + self.dimensions[0] + 1
        ymax = xmax + self.dimensions[1] + 1
        zmax = ymax + self.dimensions[2] + 1
        min_index = (xmin, xmax, ymax)
        max_index = (xmax, ymax, zmax)
        self.bounds = [[values[i] for i in range(min_index[j], max_index[j])] for j in range(0, 3)]

        dose_min_index = zmax
        dose_max_index = dose_min_index + self.nvox()
        dose = np.array(values[dose_min_index:dose_max_index])
        dose = dose.reshape(self.dimensions, order='F')

        fractional_uncertainty = np.array(values[dose_max_index:])
        fractional_uncertainty = fractional_uncertainty.reshape(self.dimensions, order='F')

        return dose, fractional_uncertainty

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
        file.write('\n'.join(' '.join(format(b, '.4f') for b in d) for d in self.bounds) + '\n')
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

    def get_max_dose_and_ijk(self) -> (float, float, tuple):
        max_index = np.unravel_index(np.argmax(self.dose), tuple(self.dimensions))

        return np.max(self.dose), self.fract_unc[max_index], max_index

    def get_min_dose(self) -> (float, float):
        """Return the minimum dose of the distribution and its fractional uncertainty

        Returns
            (float, float): the minimum dose and its fractional uncertainty
        """
        min_index = np.unravel_index(np.argmin(self.dose), tuple(self.dimensions))

        return np.min(self.dose), self.fract_unc[min_index]

    def get_min_dose_and_ijk(self) -> (float, float, tuple):
        min_index = np.unravel_index(np.argmin(self.dose), tuple(self.dimensions))

        return np.max(self.dose), self.fract_unc[min_index], min_index

    def get_point_dose_from_xyz(self, pos: Tuple[float, float, float]) -> (float, float):
        """Return the dose and its fractional uncertainty in the voxel containing position (x, y, z)

        Args:
            pos (tuple(float)): the position (x, y, z)

        Returns:
            (float, float): the dose and its fractional uncertainty
        """
        index = voxelnav.get_ijk_from_xyz(pos, self.bounds)

        return self.dose[index], self.fract_unc[index]

    def get_point_dose_from_ijk(self, index: Tuple[int, int, int]) -> (float, float):

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


def write_3ddose_to_dicom(the_dose: DoseDistribution, dicom_file=None, flip_zaxis=True):
    if dicom_file:
        dicom_dataset = pydicom.read_file(dicom_file)
        dicom_dataset.SeriesInstanceUID += "{:03d}".format(random.randint(1, 1000))
    else:
        path = join(join(EGS_TOOLS_HOME, "templates"), EMPTY_DICOM_TEMPLATE)
        dicom_dataset = pydicom.read_file(path)
        sop, study, series, frame = dicom.uid.generate_new_uids('rtdose')
        dicom_dataset.SOPInstanceUID = sop
        dicom_dataset.StudyInstanceUID = study
        dicom_dataset.SeriesInstanceUID = series
        dicom_dataset.FrameOfReferenceUID = frame

    dicom_dataset.NumberOfFrames = str(the_dose.dimensions[2])
    dicom_dataset.Rows = the_dose.dimensions[1]
    dicom_dataset.Columns = the_dose.dimensions[0]

    dx, dy, dz = voxelnav.get_voxel_size_from_ijk((0, 0, 0), the_dose.bounds)
    dicom_dataset.PixelSpacing = [str(dx * 10), str(dy * 10)]

    if flip_zaxis:
        dicom_dataset.ImagePositionPatient = [str(i * 10) for i in voxelnav.get_voxel_center_from_ijk((0, 0, -2),
                                                                                                      the_dose.bounds)]
        dicom_dataset.GridFrameOffsetVector = [str(i * -dz * 10) for i in range(0, the_dose.dimensions[2])]
    else:
        dicom_dataset.ImagePositionPatient = [str(i * 10) for i in voxelnav.get_voxel_center_from_ijk((0, 0, 0),
                                                                                                      the_dose.bounds)]
        dicom_dataset.GridFrameOffsetVector = [str(i * dz * 10) for i in range(0, the_dose.dimensions[2])]

    dose_max = the_dose.dose.max()
    num_bytes = 4
    float2int = vmap.FloatingPointToIntegerMapping(0., dose_max, num_bytes, mode='bytes')

    dicom_dataset.DoseGridScaling = str(float2int.reverse_mapping_factor)
    dose_array_in_int = float2int.float_to_integer(the_dose.dose).astype('uint32')
    dose_array = np.swapaxes(dose_array_in_int, 0, 2)  # in dicom, arrays are stored as [z][y][x] with x fastest moving
    if flip_zaxis:
        dose_array = np.flip(dose_array, 0)
    dicom_dataset.PixelData = dose_array.tostring()

    return dicom_dataset
