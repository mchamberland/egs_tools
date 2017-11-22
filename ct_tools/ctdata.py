"""ctdata

This module exports classes to handle CT image data.

class:

CTframe -- holds a single x-y CT slice as a flattened list
CTdata -- a CTdata object
"""
import numpy as np
from typing import List


def get_list_of_ctframes(ct_files):
    ctframes = []
    print('Assembling list of CT frames...')
    for (i, file) in enumerate(ct_files):
        print('Frame {0:d} out of {1:d}...'.format(i+1, len(ct_files)))
        ctframes.append(CTframe(file))
    print('Done!\n')
    return ctframes


class CTframe:
    def __init__(self, ct_file):
        self.number_of_rows = ct_file.Rows
        self.number_of_columns = ct_file.Columns
        self.row_spacing_in_cm = float(ct_file.PixelSpacing[0]) / 10
        self.column_spacing_in_cm = float(ct_file.PixelSpacing[1]) / 10
        self.image_position_in_cm = [float(x) / 10 for x in ct_file.ImagePositionPatient]
        self.rescale_intercept = float(ct_file.RescaleIntercept)
        self.rescale_slope = float(ct_file.RescaleSlope)

        self.flat_raw_values = []
        for y in range(self.number_of_rows):
            for x in range(self.number_of_columns):
                self.flat_raw_values.append(ct_file.pixel_array[y][x])


class CTdata:
    def __init__(self, ctframes: List[CTframe]):
        self.dimensions = [ctframes[0].number_of_columns, ctframes[0].number_of_rows, len(ctframes)]
        self.bounds = self.get_bounds(ctframes)
        self.image = np.array([])
        self.get_ct_data(ctframes)

    @staticmethod
    def get_bounds(ctframes: List[CTframe]):
        xbounds = []
        ybounds = []
        zbounds = []

        delta_z = ctframes[1].image_position_in_cm[2] - ctframes[0].image_position_in_cm[2]

        for frame in ctframes:
            delta_x = frame.column_spacing_in_cm
            delta_y = frame.row_spacing_in_cm

            xbounds.append(frame.image_position_in_cm[0] - delta_x / 2)
            ybounds.append(frame.image_position_in_cm[1] - delta_y / 2)
            zbounds.append(frame.image_position_in_cm[2] - delta_z / 2)

            for i in range(frame.number_of_columns):
                xbounds.append(xbounds[i] + delta_x)

            for j in range(frame.number_of_rows):
                ybounds.append(ybounds[j] + delta_y)

            for k in range(len(ctframes)):
                zbounds.append(zbounds[k] + delta_z)

        return [xbounds, ybounds, zbounds]

    def get_ct_data(self, ctframes: List[CTframe]):
        for frame in ctframes:
            temp = np.array(frame.flat_raw_values) * frame.rescale_slope + frame.rescale_intercept
            self.image = np.append(self.image, temp)
        self.image = self.image.reshape(self.dimensions, order='F')

    def nvox(self) -> int:
        """Return the total number of voxels"""
        nvox = 1
        for dims in self.dimensions:
            nvox *= dims
        return nvox
