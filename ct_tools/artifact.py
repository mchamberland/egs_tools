import voxelnav
from typing import List
from ct_tools.ctdata import CTdata


class SimpleThresholdReplacement:
    def __init__(self, threshold=200, replacement=25, low_threshold=-100, low_replacement=-50,
                 xy_search_in_voxels=10, z_search_in_voxels=2):

        self.xy_search = xy_search_in_voxels
        self.z_search = z_search_in_voxels
        self.threshold = threshold
        self.replacement = replacement
        self.low_threshold = low_threshold
        self.low_replacement = low_replacement

    def apply_str_to_seed_locations(self, ctdata: CTdata, seed_locations: List[List[float]]):
        for (i, location) in enumerate(seed_locations):
            i, j, k = voxelnav.get_ijk_from_xyz(location, ctdata.bounds)
            self.check_indices_for_out_of_bounds_search(ctdata, i, j, k)
            for delta_z in range(-self.z_search, self.z_search + 1):
                for delta_y in range(-self.xy_search, self.xy_search):
                    for delta_x in range(-self.xy_search, self.xy_search):
                        if ctdata.image[(i+delta_x, j+delta_y, k+delta_z)] > self.threshold:
                            ctdata.image[(i + delta_x, j + delta_y, k + delta_z)] = self.replacement
                        else:
                            pass

    def check_indices_for_out_of_bounds_search(self, ctdata: CTdata, i, j, k):
        if i < self.xy_search or i > ctdata.dimensions[0] - self.xy_search or \
           j < self.xy_search or j > ctdata.dimensions[1] - self.xy_search or \
           k < self.z_search or k > ctdata.dimensions[2] - self.z_search:
            raise Exception("Search around STR locations will go out of bounds!")
        else:
            pass
