import os
import numpy as np
from os.path import join


class CTConversionToTissue:
    def __init__(self, filename=None, directory="."):
        self.directory = directory
        self.media_list = []
        self.ctnum_bins = []
        if filename:
            self.read_ctconv_file(filename)

    def read_ctconv_file(self, filename='default'):
        if not filename.endswith('.ctconv'):
            filename += '.ctconv'
        path = join(self.directory, filename)
        if not os.path.exists(path):
            raise FileNotFoundError

        with open(path, 'r') as file:
            lines = file.readlines()
            lines = [line.strip() for line in lines if line.strip()]

        previous_max_ctnum = -999999
        for line in lines:
            name, density, min_ctnum, max_ctnum = line.split()
            if int(min_ctnum) < previous_max_ctnum:
                raise Exception('Min HU number for {medium} is smaller than max HU number '
                                'of previous medium.'.format(medium=name))
            if int(max_ctnum) < int(min_ctnum):
                raise Exception('Max HU number for {medium} is smaller than min HU number.'.format(medium=name))
            previous_max_ctnum = int(max_ctnum)
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

    def get_media_name_list(self):
        media_name = []
        for medium in self.media_list:
            media_name.append(medium.name)
        return media_name


class MediumInfo:
    def __init__(self, name='WATER_0.998', density=0.998, min_ctnum=-9999, max_ctnum=9999):
        self.name = name
        self.density = density
        self.min_ctnum = min_ctnum
        self.max_ctnum = max_ctnum
