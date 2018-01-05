import os
import voxelnav
import numpy as np
import egsphant.manip as egsphantmanip
from numba import jit
from os.path import join
from collections import defaultdict
from ct_tools.hu2rho import HU2rho
from egsphant.egsphant import EGSphant
from ct_tools.ct_to_tissue import CTConversionToTissue


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
        ctdata_dict = setup_ctdata_dictionary(ctdata)
        contour_path_dict = setup_contour_path_dictionary(ctdata, self.contour_info_dictionary)

        print("Converting CT data to egsphant. This may take a while...")
        # easy case: no contours, use REMAINDER
        if not self.contour_info_dictionary:
            contour = 'REMAINDER'
            loop_counter = 0
            print_counter = 10
            n = int(ctdata.nvox() / 10)
            for (index, ctnum) in np.ndenumerate(ctdata.image):
                loop_counter += 1
                if loop_counter % n == 0:
                    print("{:d}%...".format(print_counter))
                    print_counter += 10
                medium = self.tissue_converter[contour].get_medium_name_from_ctnum(ctnum)
                medium_key = egsphant.get_medium_key(medium)
                egsphant.phantom[index] = medium_key

                if self.density_instruction[contour] == 'CT':
                    egsphant.density[index] = self.density_converter.get_density_from_hu(ctnum, extrapolate)
                elif self.density_instruction[contour] == 'NOMINAL':
                    medium_index = self.tissue_converter[contour].get_medium_index_from_ctnum(ctnum)
                    density = self.tissue_converter[contour].get_medium_density(medium_index)
                    egsphant.density[index] = density
                else:
                    egsphant.density[index] = float(self.density_instruction[contour])
        else:
            egsphant = self._convert_ct_using_contours(egsphant, ctdata, ctdata_dict, contour_path_dict, extrapolate)

        print("Conversion completed! (Whew!)")
        return egsphant

    def _convert_ct_using_contours(self, egsphant, ctdata, ctdata_dict, contour_path_dict, extrapolate):
        loop_counter = 0
        print_counter = 10
        n = int(ctdata.nvox() / 10)
        for (index, value) in ctdata_dict.items():
            loop_counter += 1
            if loop_counter % n == 0:
                print("{:d}%...".format(print_counter))
                print_counter += 10

            ctnum, (x, y), k = value
            in_structure = None
            found_structure = False
            for name in self.contour_order[0:-1]:  # last contour is 'REMAINDER'
                if k in contour_path_dict[name]:
                    # TODO check if we should even bother checking the point
                    for path in contour_path_dict[name][k]:
                        if path.contains_point((x, y)):
                            in_structure = name
                            found_structure = True
                            break
                    if found_structure:
                        break

            if not found_structure:
                in_structure = 'REMAINDER'

            egsphant.phantom[index] = egsphant.inverse_key_mapping[
                self.tissue_converter[in_structure].get_medium_name_from_ctnum(ctnum)]

            if self.density_instruction[in_structure] == 'CT':
                egsphant.density[index] = self.density_converter.get_density_from_hu(ctnum, extrapolate)
            elif self.density_instruction[in_structure] == 'NOMINAL':
                egsphant.density[index] = self.tissue_converter[in_structure].get_medium_density(
                    self.tissue_converter[in_structure].get_medium_index_from_ctnum(ctnum))
            else:
                egsphant.density[index] = float(self.density_instruction[in_structure])
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


def setup_ctdata_dictionary(ctdata):
    ctdata_dict = {}
    loop_counter = 0
    print_counter = 10
    n = int(ctdata.nvox() / 10)
    xybounds = ctdata.bounds[0:2]
    print("Setting up the CT data for conversion...")
    # TODO do not calculate the pixel centre from ij for EVERY slice! They're the same!
    for (index, ctnum) in np.ndenumerate(ctdata.image):
        i, j, k = index
        loop_counter += 1
        if loop_counter % n == 0:
            print("{:d}%...".format(print_counter))
            print_counter += 10
        ctdata_dict[index] = (ctnum, voxelnav.get_pixel_center_from_ij((i, j), xybounds), k)
    print("CT data ready!")
    return ctdata_dict


def setup_contour_path_dictionary(ctdata, contour_info_dict):
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

    print("Contour data ready!")
    return contour_path_dict
