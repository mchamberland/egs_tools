import voxelnav
import numpy as np
import dicom.reader as bdr
from typing import List
from scipy.spatial import cKDTree
from matplotlib.path import Path


def get_contours_from_dicom(directory='.'):
    filenames, rt_structs = bdr.read_structure_files_in_directory(directory)
    if not rt_structs:
        return None
    if len(rt_structs) > 1:
        raise Exception("There are more than 1 RT Structure files in the directory."
                        "Please select a directory with a single RT Structure file")
    struct = rt_structs[0]
    contour_dict = {}
    for ROI in struct.StructureSetROISequence:
        label = ROI.ROIName.replace(' ', '_').lower()
        number = ROI.ROINumber
        contour_dict[number] = Contour(number=number, name=label)

    for ROI in struct.RTROIObservationsSequence:
        number = ROI.ReferencedROINumber
        if number in contour_dict:
            contour_dict[number].type = ROI.RTROIInterpretedType

    for ROI in struct.ROIContourSequence:
        number = ROI.ReferencedROINumber

        if number in contour_dict:
            if 'ROIDisplayColor' in ROI:
                contour_dict[number].colour = [float(v) for v in ROI.ROIDisplayColor]

            temp_data = []
            if 'ContourSequence' in ROI:
                for zslice in ROI.ContourSequence:
                    points_as_float_in_cm = [float(v) / 10 for v in zslice.ContourData]
                    it = iter(points_as_float_in_cm)
                    temp_data.append(list(zip(it, it, it)))

            temp_data.sort(key=lambda temp: temp[0][2])  # sort by z-slice

            zslices = []
            contour_data = []
            contour_as_path = {}
            for zslice in temp_data:
                x, y, z = list(zip(*zslice))
                contour_data.append(list(zip(x, y)))
                # some structures have multiple contours on the same slice (e.g., ribs), so first check that there are
                # no contours currently stored for this slice; otherwise, append the contour to the list
                z_value = round(z[0], 4)
                if z_value not in contour_as_path:
                    contour_as_path[z_value] = [Path(list(zip(x, y)), closed=True)]
                else:
                    contour_as_path[z_value].append(Path(list(zip(x, y)), closed=True))

                zslices.append(z[0])
            contour_dict[number].zslices = zslices
            contour_dict[number].contour_data = contour_data
            contour_dict[number].contour_as_path = contour_as_path

    contour_dict_by_label = {}
    for (key, value) in contour_dict.items():
        if contour_dict[key].type == 'BRACHY_CHANNEL':
            pass
        else:
            contour_dict_by_label[contour_dict[key].name] = contour_dict[key]

    return contour_dict_by_label


def interpolate_contours_from_dict_along_z(contour_dict, ctdata):
    slice_thickness = ctdata.voxel_size_in_cm()[2]
    interpolated_contours_dict = {}
    for label, contour in contour_dict.items():
        if contour.is_missing_slices(slice_thickness):
            interpolated_contours_dict[label] = contour.interpolate_missing_slices()
    return interpolated_contours_dict


class Contour:
    def __init__(self, contour_data=None, zslices=None, name=None, number=None, roitype=None, colour=(255, 0, 0)):
        self.name = name
        self.number = number
        self.type = roitype
        self.colour = colour
        self.zslices = zslices  # the z-slices corresponding to the contour_data, in cm
        self.contour_data = contour_data  # 2-D contour data for each zslice stored as a list of (x,y) tuples, in cm
        self.contour_as_path = None
        self.missing_slices = None
        self.is_interpolated = False
        self.pixel_indices = None
        self.zindices = None
        self.contour_mask = None

    def is_missing_slices(self, slice_thickness):
        slice_differences = np.diff(self.zslices) / slice_thickness
        missing_slices = np.around(slice_differences - 1, decimals=8)
        total_number_of_missing_slices = int(missing_slices.sum())
        if total_number_of_missing_slices > 0:
            print("Missing {} slice(s) in contour {}:".format(total_number_of_missing_slices, self.name))
            self.missing_slices = []
            for index in np.where(missing_slices > 0)[0]:
                print("Missing {} slice(s) between z = {} and z = {}".format(int(missing_slices[index]),
                                                                             self.zslices[index],
                                                                             self.zslices[index + 1]))
                self.missing_slices.append(MissingSlices(index, int(missing_slices[index]), slice_thickness))
            print("")
            return True
        else:
            return False

    def interpolate_missing_slices(self):
        if not self.missing_slices:
            print("Contours are not missing slices. Nothing to do here!")
        else:
            print("Interpolating missing slices...")
            for (i, missing_slice_group) in enumerate(self.missing_slices):
                start = missing_slice_group.start_index
                start_slice = np.array(self.contour_data[start])
                end_slice = np.array(self.contour_data[start + 1])

                # find nearest neighbours of the end slice in the first slice
                distances, indices = cKDTree(start_slice).query(end_slice)
                point_pairs = list(zip(start_slice[indices], end_slice))
                zstart = self.zslices[start]
                zend = self.zslices[start + 1]

                number_missing = missing_slice_group.number_of_missing_slices

                for n in range(number_missing):
                    zslice = self.zslices[start] + ((n + 1) * missing_slice_group.thickness)
                    interpolated_xy = []
                    for (start_xy, end_xy) in point_pairs:
                        x1, y1 = start_xy
                        x2, y2 = end_xy
                        x = (x1 * (zend - zslice) + x2 * (zslice - zstart)) / (zend - zstart)
                        y = (y1 * (zend - zslice) + y2 * (zslice - zstart)) / (zend - zstart)
                        interpolated_xy.append((x, y))
                    self.zslices = np.append(self.zslices, zslice)
                    self.contour_data.append(interpolated_xy)

            print("...Done! Sorting the newly added slices...")
            temp = list(zip(self.contour_data, self.zslices))
            temp.sort(key=lambda f: f[1])
            self.contour_data, self.zslices = zip(*temp)
            self.update_contour_as_path()
            self.missing_slices = None
            self.is_interpolated = True
            print("...Sorted!")

            return self

    def update_contour_as_path(self):
        self.contour_as_path = {}
        for index, zslice in enumerate(self.zslices):
            z_value = round(zslice, 4)
            if z_value not in self.contour_as_path:
                self.contour_as_path[z_value] = [Path(self.contour_data[index], closed=True)]
            else:
                self.contour_as_path[z_value].append(Path(self.contour_data[index], closed=True))

    def calculate_voxel_indices_from_contour_data(self, bounds3d: List[List[float]]):
        zindices = []
        pixel_indices = []
        for (index, zslice) in enumerate(self.contour_data):
            indices = []
            for (x, y) in zslice:
                i, j = voxelnav.get_ij_from_xy((x, y), bounds3d[0:2])
                indices.append((i, j))
            k = voxelnav.get_index_from_position(self.zslices[index], bounds3d[2])
            zindices.append(k)
            pixel_indices.append(indices)
        self.zindices = zindices
        self.pixel_indices = pixel_indices


class MissingSlices:
    def __init__(self, start, number, thickness):
        self.start_index = start
        self.number_of_missing_slices = number
        self.thickness = thickness


def get_extents_of_contour(zslices, paths) -> (float, float, float, float, float, float):

    zmin = min(zslices)
    zmax = max(zslices)

    xmin, ymin = 1e5, 1e5
    xmax, ymax = -1e5, -1e5
    for zslice, contours_on_slice in paths.items():
        for contour in contours_on_slice:
            xmin_tmp, ymin_tmp, xmax_tmp, ymax_tmp = contour.get_extents().extents
            if xmin_tmp < xmin:
                xmin = xmin_tmp
            if ymin_tmp < ymin:
                ymin = ymin_tmp
            if xmax_tmp > xmax:
                xmax = xmax_tmp
            if ymax_tmp > ymax:
                ymax = ymax_tmp

    return xmin, xmax, ymin, ymax, zmin, zmax
