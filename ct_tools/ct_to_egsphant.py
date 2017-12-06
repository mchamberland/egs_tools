import os
import numpy as np
import egsphant.manip as egsphantmanip
from os.path import join
from ct_tools.hu2rho import HU2rho
from egsphant.egsphant import EGSphant
from ct_tools.ct_to_tissue import CTConversionToTissue


class CTConversionToEGSphant:
    def __init__(self, filename=None, contours=None, directory='.', is_verbose=False):
        self.directory = directory
        self.contour_dictionary = contours  # instead of reading .struct files, can simply pass contours directly
        self.density_converter = None
        self.tissue_converter = {}
        self.density_instruction = {}
        self.media_list = []
        self.is_verbose = is_verbose
        if filename:
            self.read_ctscheme_file(filename)

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
            if self.contour_dictionary or contour == 'REMAINDER':
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
                # TODO read .struct file
                pass

        temp_media_list = []
        for contour, tissue_converter in self.tissue_converter.items():
            temp_media_list += tissue_converter.get_media_name_list()

        self.media_list = set(temp_media_list)

    def convert_to_egsphant(self, ctdata, extrapolate=False):
        egsphant = self.setup_egsphant(ctdata)

        print("Converting CT data to egsphant. This may take a while...")
        # easy case: no contours, use REMAINDER
        if not self.contour_dictionary:
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
            # TODO other TAS schemes
            pass

        print("Conversion completed!")
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


def is_voxel_within_max_and_min_bounds_of_contours(index, contour):
    is_within = [contour.min_indices[d] < i < contour.max_indices[d] for d, i in enumerate(index)]
    return all(is_within)