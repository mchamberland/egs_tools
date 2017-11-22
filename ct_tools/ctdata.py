"""ctdata

This module exports classes to handle CT image data.

class:

CTframe -- holds a single x-y CT slice, in raw pixel numbers
CTdata -- a CTdata object
"""
import numpy as np
from typing import List


class CTframe:
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
    def __init__(self, ctframes: List[CTframe]):
        self.dimensions = [ctframes[0].number_of_columns, ctframes[0].number_of_rows, len(ctframes)]
        self.bounds = self.get_bounds(ctframes)
        self.image = np.array([])
        self.get_ct_data(ctframes)

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

    def write_to_file(self, filename="image.ctdata"):
        with open(filename, 'wt') as ctdata:
            ctdata.write('{}\n'.format(self.number_of_media))
            for m in self.media:
                ctdata.write(m + '\n')

            dummy = ''
            for i in range(0, self.number_of_media):
                dummy += '{}\t'.format(0.25)
            dummy += '\n'
            ctdata.write(dummy)

            ctdata.write('\t'.join(str(n) for n in self.dimensions) + '\n')

            ctdata.write('\n'.join(' '.join(str(b) for b in d) for d in self.bounds))

            phantom_str = '\n'
            for k in range(0, self.dimensions[2]):
                for j in range(0, self.dimensions[1]):
                    line = ''
                    for i in range(0, self.dimensions[0]):
                        line += self.ctdata[(i, j, k)]
                    phantom_str += line + '\n'
                phantom_str += '\n'
            ctdata.write(phantom_str)



def get_list_of_ctframes(ct_files) -> List[CTframe]:
    ctframes = []
    print('Assembling list of CT frames...')
    for (i, file) in enumerate(ct_files):
        print('Frame {0:d} out of {1:d}...'.format(i+1, len(ct_files)))
        ctframes.append(CTframe(file))
    print('Done!\n')
    return ctframes
