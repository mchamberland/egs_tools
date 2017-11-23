import numpy as np
import ct_tools.ctdata as ctt
from ct_tools.ctdata import CTdata


class MetalArtifactReduction:
    def __init__(self, ctdata: CTdata, threshold=250, replacement=25, xy_search_in_voxels=10, z_search_in_voxels=2):
        self.xy_search = max(ctdata.voxel_size_in_cm()[0], ctdata.voxel_size_in_cm()[1]) * xy_search_in_voxels
        self.z_search = ctdata.voxel_size_in_cm()[2] * z_search_in_voxels
        self.threshold = threshold
        self.replacement = replacement

