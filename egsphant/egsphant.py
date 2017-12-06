"""egsphant

This module exports a single class that handles egsphant files.

class:

EGSphant -- an egsphant object
"""
import os
import gzip
import numpy as np
from typing import List


class EGSphant:
    """An egsphant object

    The class stores an egsphant object. It can read .egsphant files and the expected format is as follows:

        1) An integer representing the number of media.

        2) A list of media names (one per line).

        3) A dummy line containing ESTEPE values which are not used anymore (1 for each medium).

        4) 3 integers representing the number of voxels in each dimension.

        5) The bounds of the voxels in each dimension. X is listed first, then Y, then Z. Historically, it seems
        as if three bounds are listed on each line with a line break between the X, Y, and Z bounds.

        6) X-Y arrays of medium keys for all Z slices of the egsphant. Each medium in the list above is assigned
        a key, in order, from the dictionary of keys: '123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ' (i.e. the first medium
        in the list has key '1', the second '2', etc.)

        7) The mass density of every single voxel in the egsphant, in g / cm^3.



    Attributes:
        number_of_media (int): number of media in the egsphant
        media (list[str]): names of the media
        medium_keys (dict[str:str]): dictionary mapping the medium key (e.g.'1', '5', 'B', etc.) to the medium name
        dimensions (list[int]): number of voxels in each dimension of the distribution
        bounds (list[list[float]]): list of list of the bounds of each voxel in each dimension
        phantom (numpy.ndarray): multi-dimensional array storing the medium assigned to each voxel in the phantom
        density (numpy.ndarray): multi-dimensional array storing the mass density assigned to each voxel

    Methods:
        __init__(filename)
            Create an EGSphant instance


        read_egsphant(filename)
            Read an .egsphant file specified by its filename

        write_egsphant(filename)
            Write an egsphant to a .egsphant file specified by filename

        create_water_egsphant(bounds)
            Create a dummy water egsphant using the bounds provided


        nvox()
            Return the total number of voxels

        get_medium_key(medium)
            Return the medium key (e.g. '1', 'B', etc.) of the medium


        get_idx_from_bounds(pos, bounds)
            Get 1D index of voxel containing the position

        get_idx_from_bounds_for_region(low_pos, high_pos, bounds)
            Return a tuple of the 1D indices defining the region specified by the positions

        get_idx_of_voxels_with_medium(medium)
            Return a list of the 1-D index of the voxels containing the medium

        get_ijk_from_idx(idx):
            Return the i,j,k index of the 1-D voxel index

        get_voxel_center_from_idx(idx):
            Return the coordinates of the center of the voxel specified with 1-D index

        get_ijk_from_xyz(x, y, z):
            Return the (i, j, k) index of the voxel containing position (x, y, z)

        get_idx_from_ijk(i, j, k):
            Return the 1-D index corresponding to the (i,j,k) index

    """

    MEDIUM_KEY_STRING = '123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ'

    def __init__(self, filename=None):
        """Create an EGSphant instance

        Args:
            filename (str): filename of the egsphant file to be read

        """
        self.number_of_media = 0
        self.media = []
        self.medium_keys = dict()
        self.dimensions = []
        self.bounds = [[], [], []]
        self.phantom = np.array([0])
        self.density = np.array([0])

        if filename is not None:
            self.read_egsphant(filename)

    def read_egsphant(self, filename):
        """Read an egsphant file specified by its filename

        Args:
            filename (str): filename of the egsphant file to be read

        """
        if not os.path.exists(filename):
            return -1

        if filename.endswith('gz'):
            with gzip.open(filename, 'r') as file:
                lines = file.readlines()
                lines = [line.strip() for line in lines if line.strip()]
        else:
            with open(filename, 'r') as file:
                lines = file.readlines()
                lines = [line.strip() for line in lines if line.strip()]

        self.number_of_media = int(lines[0])
        self.media = [medium.strip() for medium in lines[1:self.number_of_media + 1]]
        self.medium_keys = dict(zip(self.MEDIUM_KEY_STRING[:len(self.media)], self.media))
        self.dimensions = list(map(int, lines[self.number_of_media+2].split()))
        del lines[0:self.number_of_media + 3]

        idx = 0
        for i in range(3):
            while len(self.bounds[i]) < self.dimensions[i] + 1:
                tmp = [float(bound) for bound in lines[idx].split()]
                self.bounds[i].extend(tmp)
                idx += 1
        del lines[0:idx]

        num_rows = self.dimensions[1] * self.dimensions[2]
        rows = [list(line) for line in lines[:num_rows]]
        phantom = [medium for row in rows for medium in row]
        self.phantom = np.array(phantom)
        self.phantom = self.phantom.reshape(self.dimensions, order='F')

        del lines[0:num_rows]

        rows = [list(map(float, line.split())) for line in lines]
        densities = [density for row in rows for density in row]
        self.density = np.array(densities)
        self.density = self.density.reshape(self.dimensions, order='F')

    def write_egsphant(self, filename):
        """Write an egsphant file specified by its filename

        Args:
            filename (str): filename of the egsphant file to be written
        """
        with open(filename, 'wt') as egsphant:
            egsphant.write('{}\n'.format(self.number_of_media))
            for m in self.media:
                egsphant.write(m + '\n')

            dummy = ''
            for i in range(0, self.number_of_media):
                dummy += '{}\t'.format(0.25)
            dummy += '\n'
            egsphant.write(dummy)

            egsphant.write('\t'.join(str(n) for n in self.dimensions) + '\n')

            egsphant.write('\n'.join(' '.join(format(b, '.4f') for b in d) for d in self.bounds))

            phantom_str = '\n'
            for k in range(0, self.dimensions[2]):
                for j in range(0, self.dimensions[1]):
                    line = ''
                    for i in range(0, self.dimensions[0]):
                        line += self.phantom[(i, j, k)]
                    phantom_str += line + '\n'
                phantom_str += '\n'
            egsphant.write(phantom_str)

            egsphant.write(' '.join('{:.8g}'.format(d) for d in self.density.flatten(order='F')) + '\n')

    @classmethod
    def create_water_egsphant(cls, bounds: List[List[float]]):
        """Create a water egsphant using the bounds provided

        Args:
            bounds (List(list[float])): the x, y, and z bounds of the egsphant to be created

        Returns
            (ESGphant): a new EGSphant instance, made of water

        """
        egsphant = cls()
        egsphant.number_of_media = 1
        egsphant.media = ['WATER_0.998']
        egsphant.medium_keys = dict([('1', 'WATER_0.998')])
        egsphant.bounds[0], egsphant.bounds[1], egsphant.bounds[2] = bounds
        egsphant.dimensions.append(len(egsphant.bounds[0]) - 1)
        egsphant.dimensions.append(len(egsphant.bounds[1]) - 1)
        egsphant.dimensions.append(len(egsphant.bounds[2]) - 1)
        egsphant.phantom = np.full(egsphant.dimensions, '1', order='F')
        egsphant.density = np.full(egsphant.dimensions, 0.998, order='F')

        return egsphant

    def nvox(self):
        """Return the total number of voxels"""
        nvox = 1
        for dims in self.dimensions:
            nvox *= dims
        return nvox

    def get_medium_key(self, medium: str):
        """Return the medium key of the medium

        Args:
            medium (str): medium name

        Returns:
            str: medium key
        """
        inverse_key_mapping = dict((v, k) for k, v in self.medium_keys.items())
        return inverse_key_mapping[medium]

    def get_index_of_voxels_with_medium(self, medium: str) -> List[int]:
        """Return a list of the 1-D index of the voxels containing 'medium'

        Args:
            medium (str): medium name

        Returns:
            list[int]: 1-D indices of voxels containing 'medium'
        """
        voxel_indices = np.where(self.phantom.flatten(order='F') == self.get_medium_key(medium))
        return voxel_indices[0].tolist()


class EGSphantFromCT(EGSphant):
    def __init__(self):
        EGSphant.__init__(self)
        self.temp_medium = Medium()
        self.voxel_assignments = []
        self.densities = []


class Medium:
    def __init__(self, name=None, density=None):
        self.name = name
        self.density = density
