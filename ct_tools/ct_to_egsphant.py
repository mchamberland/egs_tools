import os
import voxelnav
import numpy as np
import pandas as pd
import egsphant.manip as egsphantmanip
from os.path import join
from collections import defaultdict
from ct_tools.hu2rho import HU2rho
from egsphant.egsphant import EGSphant
from ct_tools.ct_to_tissue import CTConversionToTissue


# TODO delete unused contours from contour_info_dictionary so that create_contour_masks does not rely on contour_order
# TODO spit create_contour_masks out of the CTConversionToEGSphant class so it can be used by itself
class CTConversionToEGSphant:
    def __init__(self, ctscheme_filename=None, contour_dict=None, directory='.', is_verbose=False):
        self.directory = directory
        self.contour_info_dictionary = contour_dict
        self.density_converter = None
        self.tissue_converter = {}
        self.density_instruction = {}
        self.media_list = []
        self.contour_order = []
        self.is_verbose = is_verbose
        if ctscheme_filename:
            self.read_ctscheme_file(ctscheme_filename)

    def read_ctscheme_file(self, filename='default'):
        if not filename.endswith('.ctscheme'):
            filename += '.ctscheme'
        path = join(self.directory, filename)
        if not os.path.exists(path):
            return -1

        with open(path, 'r') as file:
            lines = file.readlines()
            lines = [line.strip() for line in lines if line.strip()]

        self.density_converter = HU2rho(lines[0])
        if self.is_verbose:
            print("Using {} to convert from CT numbers to density.\n".format(lines[0]))

        for line in lines[1:]:
            contour, ctconv, density_instruction = line.split()
            contour = os.path.splitext(contour)[0]
            self.contour_order.append(contour)
            if self.contour_info_dictionary or contour == 'REMAINDER':
                if self.is_verbose:
                    print("In structure {}, assign tissues using:".format(contour))
                    if not ctconv.endswith('.ctconv'):
                        print(ctconv + '.ctconv')
                    else:
                        print(ctconv)
                    if density_instruction == 'CT':
                        density_string = 'densities derived from CT numbers'
                    elif density_instruction == 'NOMINAL':
                        density_string = 'nominal densities'
                    else:
                        density_string = 'a density of {:f} g/cm^3'.format(float(density_instruction))
                    print("and {}\n".format(density_string))

                self.density_instruction[contour] = density_instruction
                self.tissue_converter[contour] = CTConversionToTissue(ctconv)
            else:
                pass

        temp_media_list = []
        for contour, tissue_converter in self.tissue_converter.items():
            temp_media_list += tissue_converter.get_media_name_list()

        self.media_list = set(temp_media_list)

    def convert_to_egsphant(self, ctdata, extrapolate=False):
        egsphant = self.setup_egsphant(ctdata)

        print("Converting CT data to egsphant...")
        if not self.contour_info_dictionary:
            mask_dict = {'REMAINDER': np.full(ctdata.dimensions, True, order='F')}
            egsphant = self._convert_using_contour_masks(egsphant, ctdata, mask_dict, extrapolate)
        else:
            contour_path_dict = setup_contour_path_dictionary(ctdata, self.contour_info_dictionary)
            print("Creating masks from contours...")
            mask_dict = self.create_contour_masks(ctdata, contour_path_dict)
            adjusted_mask_dict = self.adjust_contour_masks_by_priorities(mask_dict)
            print("Creating the egsphant...")
            egsphant = self._convert_using_contour_masks(egsphant, ctdata, adjusted_mask_dict, extrapolate)

        print("Conversion completed!")
        return egsphant

    def _convert_using_contour_masks(self, egsphant, ctdata, mask_dict, extrapolate):
        medium_key_mapping = pd.Series(egsphant.inverse_key_mapping)
        flat_phantom = egsphant.phantom.flatten(order='F')
        flat_density = egsphant.density.flatten(order='F')
        flat_image = ctdata.image.flatten(order='F')
        for contour, mask in mask_dict.items():
            flat_mask = mask.flatten(order='F')
            medium = self.tissue_converter[contour].get_medium_name_from_ctnum(flat_image[flat_mask])
            flat_phantom[flat_mask] = np.array(medium_key_mapping[medium].tolist())

            if self.density_instruction[contour] == 'CT':
                flat_density[flat_mask] = self.density_converter.get_densities_from_hu(flat_image[flat_mask],
                                                                                       extrapolate)
            elif self.density_instruction[contour] == 'NOMINAL':
                flat_density[flat_mask] = np.array(self.tissue_converter[contour].media_density_series[medium].tolist())
            else:
                flat_density[flat_mask] = float(self.density_instruction[contour])

        egsphant.phantom = flat_phantom.reshape(ctdata.dimensions, order='F')
        egsphant.density = flat_density.reshape(ctdata.dimensions, order='F')
        return egsphant

    def setup_egsphant(self, ctdata):
        egsphant = EGSphant()
        egsphant.bounds = ctdata.bounds
        egsphant.dimensions = ctdata.dimensions
        egsphant.phantom = np.empty(ctdata.dimensions, dtype=str, order='F')
        egsphant.density = np.zeros(ctdata.dimensions, order='F')
        for medium in self.media_list:
            egsphantmanip.add_medium(egsphant, medium)

        return egsphant

    def create_contour_masks(self, ctdata, contour_path_dict):
        contour_mask_dict = defaultdict(dict)
        total_cumulative_mask = np.zeros(ctdata.dimensions, dtype=bool, order='F')
        for name in self.contour_order[0:-1]:
            mask_array = np.zeros(ctdata.dimensions, dtype=bool, order='F')
            for k in range(ctdata.dimensions[2]):
                if k in contour_path_dict[name]:
                    temp_cumulative_mask = np.zeros(ctdata.dimensions[0:2], dtype=bool, order='F')
                    for path in contour_path_dict[name][k]:
                        temp_mask = path.contains_points(ctdata.pixel_centre_coordinates)
                        temp_cumulative_mask = temp_cumulative_mask | temp_mask.reshape(ctdata.dimensions[0:2],
                                                                                        order='F')
                    mask_array[:, :, k] = temp_cumulative_mask
            contour_mask_dict[name] = mask_array
            total_cumulative_mask = total_cumulative_mask | mask_array

        contour_mask_dict['REMAINDER'] = np.invert(total_cumulative_mask)

        return contour_mask_dict

    def adjust_contour_masks_by_priorities(self, contour_mask_dict):
        total_cumulative_mask = contour_mask_dict[self.contour_order[0]]
        for name in self.contour_order[1:-1]:  # nothing to do for first and last contours (last is REMAINDER)
            voxels_left = np.invert(total_cumulative_mask)
            contour_mask_dict[name] = voxels_left & contour_mask_dict[name]
            total_cumulative_mask = total_cumulative_mask | contour_mask_dict[name]

        return contour_mask_dict


def setup_contour_path_dictionary(ctdata, contour_info_dict):
    print("Preparing contour data...")
    contour_path_dict = defaultdict(dict)
    zbounds = ctdata.bounds[2]
    for name, contour in contour_info_dict.items():
        for zslice in contour.zslices:
            k = voxelnav.get_index_from_position(zslice, zbounds)
            # some structures have multiple contours on the same slice (e.g., ribs), so first check that there are no
            # contours currently stored for this slice; otherwise, append the contour to the list
            if k not in contour_path_dict:
                contour_path_dict[name][k] = [contour_info_dict[name].contour_as_path[round(zslice, 4)]]
            else:
                contour_path_dict[name][k].append(contour_info_dict[name].contour_as_path[round(zslice, 4)])

    return contour_path_dict


