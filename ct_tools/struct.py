import numpy as np
import brachy_dicom.reader as bdr


def get_contours_from_dicom(directory='.'):
    filenames, rt_structs = bdr.read_structure_files_in_directory(directory)
    if len(rt_structs) > 1:
        raise Exception("There are more than 1 RT Structure files in the directory."
                        "Please select a directory with a single RT Structure file")
    struct = rt_structs[0]
    contour_list = []
    for ROI in struct.RTROIObservationsSequence:
        if ROI.RTROIInterpretedType == "BRACHY_CHANNEL":
            continue

        label = ROI.ROIObservationLabel
        number = ROI.ReferencedROINumber
        colour = [float(v) for v in struct.ROIContourSequence[number].ROIDisplayColor]
        temp_data = []
        for zslice in struct.ROIContourSequence[number].ContourSequence:
            points_as_float_in_cm = [float(v) / 10 for v in zslice.ContourData]
            it = iter(points_as_float_in_cm)
            temp_data.append(list(zip(it, it, it)))

        temp_data.sort(key=lambda temp: temp[0][2])  # sort by z-slice
        zslices = []
        contour_data = []
        for zslice in temp_data:
            x, y, z = list(zip(*zslice))
            contour_data.append(list(zip(x, y)))
            zslices.append(z[0])

        contour_list.append(Contour(np.array(contour_data), np.array(zslices), label, number, colour))

    return contour_list


class Contour:
    def __init__(self, contour_data, zslices, name='', number=None, colour=(255, 0, 0)):
        self.name = name
        self.number = number
        self.colour = colour
        self.zslices = zslices  # the z-slices corresponding to the contour_data, in cm
        self.contour_data = contour_data  # 2-D contour data stored as (x,y) tuples, in cm
        self.missing_slices = None
        self.is_interpolated = False

    def is_missing_slices(self, slice_thickness):
        slice_differences = np.diff(self.zslices) / slice_thickness
        missing_slices = slice_differences - 1
        total_number_of_missing_slices = int(round(missing_slices.sum()))
        if total_number_of_missing_slices > 0:
            print("Missing {} slice(s) in contour {}".format(total_number_of_missing_slices, self.name))
            self.missing_slices = []
            for index in np.where(missing_slices > 0):
                self.missing_slices.append(MissingSlices(self.zslices[index][0], self.zslices[index+1][0],
                                                         int(missing_slices[index][0]), slice_thickness))
            return True
        else:
            return False


class MissingSlices:
    def __init__(self, start, end, number, thickness):
        self.start_slice = start
        self.end_slice = end
        self.number_of_missing_slices = number
        self.thickness = thickness
