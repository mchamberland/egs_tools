import os
import numpy as np
import pandas as pd
from os.path import join


class CTConversionToTissue:
    def __init__(self, filename=None, directory="."):
        self.directory = directory
        self.media_list = []
        self.media_series = pd.Series()
        self.media_density_dict = {}
        self.media_density_series = pd.Series()
        self.ctnum_bins = []
        if filename:
            self.read_ctconv_file(filename)
            self.media_series = self.build_media_series()

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
            self.media_density_dict[name] = float(density)
            self.media_list.append(medium)
            self.ctnum_bins.append(int(min_ctnum))

        self.media_density_series = pd.Series(self.media_density_dict)
        self.ctnum_bins.append(self.media_list[-1].max_ctnum)

    def get_medium_density(self, index):
        return self.media_list[index].density

    def get_medium_index_from_ctnum(self, ctnum):
        return np.searchsorted(self.ctnum_bins, ctnum) - 1

    def get_medium_name_from_ctnum(self, ctnum):
        the_indices = np.searchsorted(self.ctnum_bins, ctnum) - 1
        if np.any(the_indices >= len(self.media_list)):
            the_indices[the_indices >= len(self.media_list)] = len(self.media_list) - 1
        return np.array(self.media_series.reindex(the_indices).tolist())

    def get_media_name_list(self):
        media_name = []
        for medium in self.media_list:
            media_name.append(medium.name)
        return media_name

    def build_media_series(self):
        temp_dict = {}
        for index, medium in enumerate(self.media_list):
            temp_dict[index] = medium.name
        return pd.Series(temp_dict)


class MediumInfo:
    def __init__(self, name='WATER_0.998', density=0.998, min_ctnum=-9999, max_ctnum=9999):
        self.name = name
        self.density = density
        self.min_ctnum = min_ctnum
        self.max_ctnum = max_ctnum
