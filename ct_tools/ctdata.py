"""ctdata

This module exports classes to handle CT image data.

class:

CTframe -- holds a single x-y CT slice, in raw pixel numbers
CTdata -- a CTdata object
"""
import os
import gzip
import voxelnav
import numpy as np
import brachy_dicom.reader as bdr
from typing import List, Tuple


class CTframe:
    # note that the axes of the DICOM pixel array are switched when read because the DICOM array is in the order [y][x]
    # instead of [x][y]
    def __init__(self, ct_file):
        self.number_of_rows = ct_file.Rows
        self.number_of_columns = ct_file.Columns
        self.row_spacing_in_cm = float(ct_file.PixelSpacing[0]) / 10
        self.column_spacing_in_cm = float(ct_file.PixelSpacing[1]) / 10
        self.image_position_in_cm = [float(x) / 10 for x in ct_file.ImagePositionPatient]
        self.rescale_intercept = float(ct_file.RescaleIntercept)
        self.rescale_slope = float(ct_file.RescaleSlope)
        self.pixel_array = np.swapaxes(ct_file.pixel_array, 0, 1)


class CTdata:
    def __init__(self, ctframes: List[CTframe]=None):
        if ctframes:
            self.dimensions = [ctframes[0].number_of_columns, ctframes[0].number_of_rows, len(ctframes)]
            self.bounds = self.get_bounds(ctframes)
            self.image = np.array([])
            self.get_ct_data(ctframes)
        else:
            self.dimensions = [0, 0, 0]
            self.bounds = [[], [], []]
            self.image = np.array([])

    @staticmethod
    def get_bounds(ctframes: List[CTframe]) -> List[List[float]]:
        xbounds = []
        ybounds = []
        zbounds = []

        delta_z = ctframes[1].image_position_in_cm[2] - ctframes[0].image_position_in_cm[2]

        delta_x = ctframes[0].column_spacing_in_cm
        delta_y = ctframes[0].row_spacing_in_cm

        xbounds.append(ctframes[0].image_position_in_cm[0] - delta_x / 2)
        ybounds.append(ctframes[0].image_position_in_cm[1] - delta_y / 2)
        zbounds.append(ctframes[0].image_position_in_cm[2] - delta_z / 2)

        for i in range(ctframes[0].number_of_columns):
            xbounds.append(xbounds[i] + delta_x)

        for j in range(ctframes[0].number_of_rows):
            ybounds.append(ybounds[j] + delta_y)

        for k in range(len(ctframes)):
            zbounds.append(zbounds[k] + delta_z)

        return [xbounds, ybounds, zbounds]

    def get_ct_data(self, ctframes: List[CTframe]):
        pixel_arrays = []
        for frame in ctframes:
            pixel_arrays.append(frame.pixel_array * frame.rescale_slope + frame.rescale_intercept)
        self.image = np.stack(pixel_arrays, axis=-1)

    def nvox(self) -> int:
        """Return the total number of voxels"""
        nvox = 1
        for dims in self.dimensions:
            nvox *= dims
        return nvox

    def voxel_size_in_cm(self) -> (float, float, float):
        dx = self.bounds[0][1] - self.bounds[0][0]
        dy = self.bounds[1][1] - self.bounds[1][0]
        dz = self.bounds[2][1] - self.bounds[2][0]
        return dx, dy, dz

    def image_size_in_cm(self) -> (float, float, float):
        dx = self.bounds[0][-1] - self.bounds[0][0]
        dy = self.bounds[1][-1] - self.bounds[1][0]
        dz = self.bounds[2][-1] - self.bounds[2][0]
        return dx, dy, dz

    def get_ct_number_from_xyz(self, pos: Tuple[float, float, float]) -> float:
        index = voxelnav.get_ijk_from_xyz(pos, self.bounds)
        return self.image[index]

    def write_to_file(self, filename="image.ctdata"):
        if filename.endswith('gz'):
            file = gzip.open(filename, 'wt')
        else:
            file = open(filename, 'wt')

        file.write('\t'.join(str(n) for n in self.dimensions) + '\n')
        file.write('\n'.join(' '.join(str(b) for b in d) for d in self.bounds))

        image_str = '\n'
        for k in range(0, self.dimensions[2]):
            for j in range(0, self.dimensions[1]):
                line = ''
                for i in range(0, self.dimensions[0]):
                    line += str(self.image[(i, j, k)]) + ' '
                image_str += line + '\n'
            image_str += '\n'
        file.write(image_str)
        file.close()


def crop_ctdata_to_bounds(original: CTdata, the_bounds, include_partial_voxel=False) -> CTdata:
    if len(the_bounds) > 1:
        xi, xf, yi, yf, zi, zf = the_bounds
    else:
        xi = yi = zi = (-1) * the_bounds[0] / 2.
        xf = yf = zf = the_bounds[0] / 2.

    x_index = voxelnav.get_region_indices_from_low_and_hi_pos(xi, xf, original.bounds[0], include_partial_voxel)
    y_index = voxelnav.get_region_indices_from_low_and_hi_pos(yi, yf, original.bounds[1], include_partial_voxel)
    z_index = voxelnav.get_region_indices_from_low_and_hi_pos(zi, zf, original.bounds[2], include_partial_voxel)

    xb_slice = slice(x_index[0], x_index[1] + 2)
    yb_slice = slice(y_index[0], y_index[1] + 2)
    zb_slice = slice(z_index[0], z_index[1] + 2)

    cropped = CTdata()

    cropped.bounds[0] = original.bounds[0][xb_slice]
    cropped.bounds[1] = original.bounds[1][yb_slice]
    cropped.bounds[2] = original.bounds[2][zb_slice]

    x_slice = slice(x_index[0], x_index[1] + 1)
    y_slice = slice(y_index[0], y_index[1] + 1)
    z_slice = slice(z_index[0], z_index[1] + 1)

    cropped.image = original.image[x_slice, y_slice, z_slice]

    cropped.dimensions[0] = len(cropped.bounds[0]) - 1
    cropped.dimensions[1] = len(cropped.bounds[1]) - 1
    cropped.dimensions[2] = len(cropped.bounds[2]) - 1

    return cropped


def crop_ctdata_to_slice_indices(original: CTdata, slice1: int, slice2: int) -> CTdata:
    cropped = CTdata()
    if slice1 > slice2:
        tmp = slice1
        slice1 = slice2
        slice2 = tmp

    zb_slice = slice(slice1, slice2 + 2)
    cropped.bounds[0] = original.bounds[0]
    cropped.bounds[1] = original.bounds[1]
    cropped.bounds[2] = original.bounds[2][zb_slice]

    z_slice = slice(slice1, slice2 + 1)
    cropped.image = original.image[:, :, z_slice]

    cropped.dimensions[0] = original.dimensions[0]
    cropped.dimensions[1] = original.dimensions[1]
    cropped.dimensions[2] = len(cropped.bounds[2]) - 1

    return cropped


def read_ctdata(filename):
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

    ctdata = CTdata()
    ctdata.dimensions = list(map(int, lines[0].split()))
    del lines[0]

    idx = 0
    for i in range(3):
        while len(ctdata.bounds[i]) < ctdata.dimensions[i] + 1:
            tmp = [float(bound) for bound in lines[idx].split()]
            ctdata.bounds[i].extend(tmp)
            idx += 1
    del lines[0:idx]

    image = [float(values) for row in lines for values in row.split()]
    ctdata.image = np.array(image)
    ctdata.image = ctdata.image.reshape(ctdata.dimensions, order='F')

    del lines
    return ctdata


def get_list_of_ctframes(ct_files, verbose=False) -> List[CTframe]:
    ctframes = []
    if verbose:
        print('Assembling list of CT frames...')
    for (i, file) in enumerate(ct_files):
        if verbose:
            print('Frame {0:d} out of {1:d}...'.format(i+1, len(ct_files)))
        ctframes.append(CTframe(file))
    if verbose:
        print('Done!\n')
    return ctframes


def get_ctdata_from_dicom(directory='.') -> CTdata:
    ct_names, ct_files = bdr.get_sorted_ct_files_in_directory(directory)
    ct_frames = get_list_of_ctframes(ct_files)
    return CTdata(ct_frames)
