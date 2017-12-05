import os
import numpy as np
from os.path import join, isfile, splitext


class CTConversionToTissue:
    def __init__(self, filename=None, directory="."):
        self.directory = directory
        self.media_list = []
        self.ctnum_bins = []
        if filename:
            self.read_ctconv_file(filename)

    def read_ctconv_file(self, filename='default'):
        # TODO check CT number limits and that CT numbers are sorted in ascending order
        path = join(self.directory, filename + '.ctconv')
        if not os.path.exists(path):
            return -1

        with open(path, 'r') as file:
            lines = file.readlines()
            lines = [line.strip() for line in lines if line.strip()]

        for line in lines:
            name, density, min_ctnum, max_ctnum = line.split()
            medium = MediumInfo(name, float(density), int(min_ctnum), int(max_ctnum))
            self.media_list.append(medium)
            self.ctnum_bins.append(int(min_ctnum))

        self.ctnum_bins.append(self.media_list[-1].max_ctnum)

    def get_medium_density(self, index):
        return self.media_list[index].density

    def get_medium_index_from_ctnum(self, ctnum):
        return np.searchsorted(self.ctnum_bins, ctnum) - 1

    def get_medium_name_from_ctnum(self, ctnum):
        return self.media_list[np.searchsorted(self.ctnum_bins, ctnum) - 1].name


class MediumInfo:
    def __init__(self, name='WATER_0.998', density=0.998, min_ctnum=-9999, max_ctnum=9999):
        self.name = name
        self.density = density
        self.min_ctnum = min_ctnum
        self.max_ctnum = max_ctnum
