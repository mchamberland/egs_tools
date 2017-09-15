"""EGSinp

This module exports a single class that handles creation of egsinp files for egs_brachy.

class:

EGSinp -- egsinp file for egs_brachy
"""
import os.path


class EGSinp:
    """egsinp file for egs_brachy


    Attributes:

    Methods:
        __init__(filename)
            Create a DoseDistribution instance from a 3ddose filename


    """

    SOURCE_RADIONUCLIDE = {'OncoSeed_6711': 'I125',
                            'TheraSeed_200': 'Pd103',
                            'microSelectron-v2': 'Ir192',
                            'GammaMedPlus': 'Ir192'}

    RADIONUCLIDE_SPECTRUM = {'I125': 'I125_NCRP_line',
                             'Pd103': 'Pd103_NNDC_2.6_line',
                             'Ir192': 'Ir192_NNDC_2.6_line'}

    def __init__(self, filename="egsinp", source_model="OncoSeed_6711", path_type="relative"):
        """Create a skeleton to build an egsinp file for egs_brachy"""
        self.source_model = source_model
        self.root = ""
        if path_type == "absolute":
            self.root = os.path.expandvars("$EGS_HOME")
        self.input_file = open(filename + '.egsinp', 'w')

    def run_control(self, ncase=1e6, nbatch=1, nchunk=1, calculation='first', geometry_error_limit=2500,
                    egsdat_file_format='gzip'):
        start_delimiter = ":start run control:"
        stop_delimiter = ":stop run control:\n"
        ncase_str = "ncase = " + str(ncase)
        nbatch_str = "nbatch = " + str(nbatch)
        nchunk_str = "nchunk = " + str(nchunk)
        calc_str = "calculation = " + calculation
        geom_error_str = "geometry error limit = " + str(geometry_error_limit)
        egsdat_str = "egsdat file format = " + egsdat_file_format

        input_block = "{0}\n{1}\n{2}\n{3}\n{4}\n{5}\n{6}\n{7}\n".format(start_delimiter, ncase_str, nbatch_str,
                                                                        nchunk_str, calc_str, geom_error_str,
                                                                        egsdat_str, stop_delimiter)
        self.input_file.write(input_block)

    def run_mode(self, run_mode="normal", single_generator="yes"):
        start_delimiter = ":start run mode:"
        stop_delimiter = ":stop run mode:\n"
        mode_str = "run mode = " + run_mode
        gen_str = "single generator = " + single_generator
        input_block = "{0}\n{1}\n{2}\n{3}\n".format(start_delimiter, mode_str, gen_str, stop_delimiter)
        self.input_file.write(input_block)

    def transport_parameters(self, transport_file="low_energy_default"):
        input_block = "include file = " + self.root + "lib/transport/" + transport_file + "\n\n"
        self.input_file.write(input_block)

    def variance_reduction(self, do_recycling=True, recycling_number=1, deactivate_global_range_rejection=False,
                           global_range_rejection_max_total_energy=None):
        # TODO add source range rejection arguments
        start_delimiter = ":start variance reduction:"
        stop_delimiter = ":stop variance reduction:\n"

        if do_recycling:
            recycling_str = ":start particle recycling:\n"
            recycling_str += "times to reuse recycled particles = " + str(recycling_number) + "\n"
            recycling_str += "rotate recycled particles = yes\n"
            recycling_str += ":stop particle recycling:\n"
        else:
            recycling_str = ""

        range_rejection_str = ""
        if global_range_rejection_max_total_energy is not None:
            range_rejection_str += "global range rejection max energy = " +\
                                     str(global_range_rejection_max_total_energy) + "\n"
        if deactivate_global_range_rejection:
            range_rejection_str += "global range rejection = no\n"

        input_block = "{0}\n{1}\n{2}\n{3}\n".format(start_delimiter, recycling_str, range_rejection_str, stop_delimiter)
        self.input_file.write(input_block)

    def scoring_options(self, score_tracklength=True, score_edep=False, score_scatter=False,
                        muen_file="brachy_xcom_1.5MeV.muendat", muen_for_media="WATER_0.998", output_voxel_info=False,
                        dose_scaling_factor=None):
        start_delimiter = ":start scoring options:"
        stop_delimiter = ":stop scoring options:\n"

        if score_tracklength:
            track_str = ""
        else:
            track_str = "score tracklength dose = no"

        if score_edep:
            edep_str = "score energy deposition = yes"
        else:
            edep_str = ""

        if score_scatter:
            scatter_str = "score scatter dose = yes"
        else:
            scatter_str = ""

        muendat_str = "muen file = " + self.root + "lib/muen/" + muen_file
        muen_media_str = "muen for media = " + muen_for_media

        if output_voxel_info:
            voxel_str = "output voxel info files = yes"
        else:
            voxel_str = ""

        if dose_scaling_factor is not None:
            scaling_str = "dose scaling factor = " + str(dose_scaling_factor)
        else:
            scaling_str = ""

        input_block = "{0}\n{1}\n{2}\n{3}\n{4}\n{5}\n{6}\n{7}\n{8}\n".format(start_delimiter, track_str, edep_str,
                                                                             scatter_str, muendat_str, muen_media_str,
                                                                             voxel_str, scaling_str, stop_delimiter)
        self.input_file.write(input_block)

    def media_definition(self, ae=0.512, ue=2.012, ap=0.001, up=1.500):
        start_delimiter = ":start media definition:"
        stop_delimiter = ":stop media definition:\n"
        ae_str = "AE = " + str(ae)
        ue_str = "UE = " + str(ue)
        ap_str = "AP = " + str(ap)
        up_str = "UP = " + str(up)
        material_str = "material data file = " + self.root + "lib/media/material.dat"
        input_block = "{0}\n{1}\n{2}\n{3}\n{4}\n{5}\n".format(start_delimiter, ae_str, ue_str, ap_str, up_str,
                                                              material_str, stop_delimiter)
        self.input_file.write(input_block)

