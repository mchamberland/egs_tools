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

    SOURCE_RADIONUCLIDE = {'OncoSeed_6711': 'I125_LDR',
                           'TheraSeed_200': 'Pd103_LDR',
                           'microSelectron-v2': 'Ir192_HDR',
                           'GammaMedPlus': 'Ir192_HDR'}

    RADIONUCLIDE_SPECTRUM = {'I125_LDR': 'I125_NCRP_line',
                             'Pd103_LDR': 'Pd103_NNDC_2.6_line',
                             'Ir192_HDR': 'Ir192_NNDC_2.6_line'}

    MEAN_LIFETIME_IN_HOURS = {'I125_LDR': 2056.7,
                              'Pd103_LDR': 588.3,
                              'Ir192_HDR': 2556.3}

    # in units of Gy cm^2 hist^-1
    SOURCE_AIR_KERMA_STRENGTH_PER_HISTORY = {'OncoSeed_6711': 3.7651E-14,
                                             'TheraSeed_200': 6.4255E-14,
                                             'microSelectron-v2': 1.1517E-13,
                                             'GammaMedPlus': 1.1592E-13}

    def __init__(self, filename="egsinp", source_model="OncoSeed_6711", path_type="relative"):
        """Create a skeleton to build an egsinp file for egs_brachy"""
        self.root = ""
        self.run_mode = ""
        self.source_model = source_model
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

        input_block = "{0}\n\t{1}\n\t{2}\n\t{3}\n\t{4}\n\t{5}\n\t{6}\n{7}\n".format(start_delimiter, ncase_str,
                                                                                    nbatch_str, nchunk_str, calc_str,
                                                                                    geom_error_str, egsdat_str,
                                                                                    stop_delimiter)
        self.input_file.write(input_block)

    def run_mode(self, run_mode="normal", single_generator="yes"):
        start_delimiter = ":start run mode:"
        stop_delimiter = ":stop run mode:\n"
        mode_str = "run mode = " + run_mode
        gen_str = "single generator = " + single_generator
        input_block = "{0}\n\t{1}\n\t{2}\n{3}\n".format(start_delimiter, mode_str, gen_str, stop_delimiter)
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
            recycling_str = ":start particle recycling:\n\t"
            recycling_str += "times to reuse recycled particles = " + str(recycling_number) + "\n\t"
            recycling_str += "rotate recycled particles = yes\n"
            recycling_str += ":stop particle recycling:\n"
        else:
            recycling_str = ""

        range_rejection_str = ""
        if global_range_rejection_max_total_energy is not None:
            range_rejection_str += "global range rejection max energy = " +\
                                     str(global_range_rejection_max_total_energy) + "\n\t"
        if deactivate_global_range_rejection:
            range_rejection_str += "global range rejection = no\n"

        input_block = "{0}\n\t{1}\n\t{2}\n{3}\n".format(start_delimiter, recycling_str, range_rejection_str,
                                                        stop_delimiter)
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

        input_block = "{0}\n\t{1}\n\t{2}\n\t{3}\n\t{4}\n\t{5}\n\t{6}\n\t{7}\n{8}\n".format(start_delimiter, track_str,
                                                                                           edep_str, scatter_str,
                                                                                           muendat_str, muen_media_str,
                                                                                           voxel_str, scaling_str,
                                                                                           stop_delimiter)
        self.input_file.write(input_block)

    def media_definition(self, ae=0.512, ue=2.012, ap=0.001, up=1.500):
        start_delimiter = ":start media definition:"
        stop_delimiter = ":stop media definition:\n"
        ae_str = "AE = " + str(ae)
        ue_str = "UE = " + str(ue)
        ap_str = "AP = " + str(ap)
        up_str = "UP = " + str(up)
        material_str = "material data file = " + self.root + "lib/media/material.dat"
        input_block = "{0}\n\t{1}\n\t{2}\n\t{3}\n\t{4}\n\t{5}\n{6}\n".format(start_delimiter, ae_str, ue_str, ap_str,
                                                                             up_str, material_str, stop_delimiter)
        self.input_file.write(input_block)

    def source(self, source_type="isotropic", transformations="single_seed_at_origin", phsp=None):
        start_delimiter = ":start source definition:\n\t:start source:\n\t\t"
        stop_delimiter = ":stop source definition:\n"
        stop_source_delimiter = ":stop source:\n\t"
        source_str = "name = " + self.source_model + "\n\t\t"
        charge_str = ""
        phsp_header_str = ""
        include_str = ""

        if source_type == "isotropic":
            library_str = "library = egs_isotropic_source\n\t\t"
            charge_str = "charge = 0\n\t\t"
            include_str = "include file = " + self.root + "lib/geometry/sources/" +\
                          self.SOURCE_RADIONUCLIDE[self.source_model] + "/" +\
                          self.source_model + "/" + self.source_model + ".shape\n\t\t"
        elif phsp is not None:
            library_str = "library = eb_iaeaphsp_source\n\t\t"
            phsp_header_str = "header file = " + self.root + "lib/phsp/" + phsp + ".IAEAheader\n\t"
        else:
            library_str = ""

        simulation_source_str = "simulation source = " + self.source_model + "\n"
        start_transformations_delimiter = ":start transformations:\n\t\t"
        stop_transformations_delimiter = ":stop transformations:\n\t"
        transformations_str = "include file = " + self.root + "lib/geometry/transformations/" + transformations + "\n\t"

        start_spectrum_delimiter = ":start spectrum:\n\t\t\t"
        stop_spectrum_delimiter = ":stop spectrum:\n\t"
        spectrum_type_str = "type = tabulated spectrum\n\t\t\t"
        spectrum_file_str = "spectrum file = " + self.root + "lib/spectra/{0}.spectrum\n\t\t".format(
            self.RADIONUCLIDE_SPECTRUM[self.source_model])

        input_block = "{0}{1}{2}{3}{4}{5}{6}{7}{8}{9}{10}{11}{12}{13}{14}{15}".format(start_delimiter, source_str,
                                                                                      library_str, charge_str,
                                                                                      include_str, phsp_header_str,
                                                                                      start_spectrum_delimiter,
                                                                                      spectrum_type_str,
                                                                                      spectrum_file_str,
                                                                                      stop_spectrum_delimiter,
                                                                                      stop_source_delimiter,
                                                                                      start_transformations_delimiter,
                                                                                      transformations_str,
                                                                                      stop_transformations_delimiter,
                                                                                      simulation_source_str,
                                                                                      stop_delimiter)
        self.input_file.write(input_block)

    def volume_correction(self, correction="correct", density=1e8):
        start_delimiter = ":start volume correction:\n\t:start source volume correction:\n\t\t"
        correction_str = "correction type = {0}\n\t\t".format(correction)
        density_str = "density of random points (cm^-3) = {0}\n\t\t".format(str(density))
        boundary_str = "include file = " + self.root + "lib/geometry/sources/{0}/{1}/boundary.shape\n\t".format(
            self.SOURCE_RADIONUCLIDE[self.source_model], self.source_model)
        stop_delimiter = ":stop source volume correction:\n:stop volume correction:\n\n"
        input_block = "{0}{1}{2}{3}{4}".format(start_delimiter, correction_str, density_str, boundary_str,
                                               stop_delimiter)
        self.input_file.write(input_block)

    def geometry(self, box=None, phantom=None, egsphant=None, transformations="single_seed_at_origin",
                 discovery_action="discover", density=1e8, superposition=False):
        start_delimiter = ":start geometry definition:"

        start_geometry_delimiter = ":start geometry:"
        stop_geometry_delimiter = ":stop geometry:"

        box_block = ""
        if box is not None:
            name_str = "name = box"
            library_str = "library = egs_glib"
            include_str = "include file = " + self.root + "lib/geometry/phantoms/{0}.geom".format(box)
            box_block = "{0}\n\t\t{1}\n\t\t{2}\n\t\t{3}\n\t{4}\n".format(start_geometry_delimiter, name_str,
                                                                         library_str, include_str,
                                                                         stop_geometry_delimiter)
        phantom_block = ""
        if phantom is not None:
            name_str = "name = phantom"
            library_str = "library = egs_glib"
            include_str = "include file = " + self.root + "lib/geometry/phantoms/{0}.geom".format(phantom)
            phantom_block = "{0}\n\t\t{1}\n\t\t{2}\n\t\t{3}\n\t{4}\n".format(start_geometry_delimiter, name_str,
                                                                             library_str, include_str,
                                                                             stop_geometry_delimiter)
        egsphant_block = ""
        if egsphant is not None:
            name_str = "name = phantom"
            library_str = "library = egs_glib"
            type_str = "type = egsphant"
            egsphant_str = "egsphant file = " + self.root +\
                           "lib/geometry/phantoms/egsphant/{0}.egsphant".format(egsphant)
            egsphant_block = "{0}\n\t\t{1}\n\t\t{2}\n\t\t{3}\n\t{4}\n".format(start_geometry_delimiter, name_str,
                                                                              library_str, type_str, egsphant_str,
                                                                              stop_geometry_delimiter)

        name_str = "name = seed"
        library_str = "library = egs_glib"
        include_str = "include file = " + self.root + "lib/geometry/sources/{0}/{1}/{1}.shape\n\t".format(
            self.SOURCE_RADIONUCLIDE[self.source_model], self.source_model)
        seed_block = "{0}\n\t\t{1}\n\t\t{2}\n\t\t{3}\n\t{4}\n".format(start_geometry_delimiter, name_str,
                                                                      library_str, include_str,
                                                                      stop_geometry_delimiter)

        name_str = "name = phantom_and_seeds"
        library_str = "egs_autoenvelope"
        type_str = ""
        if superposition is True:
            type_str = "EGS_ASwitchedEnvelope"
        base_geometry_str = "base geometry = phantom"

        start_inscribed_geometry_delimiter = ":start inscribed geometry:"
        stop_inscribed_geometry_delimiter = ":stop inscribed geometry:"
        inscribed_geom_str = "inscribed geometry name = seed"

        start_transformations_delimiter = ":start transformations:"
        stop_transformations_delimiter = ":stop transformations:"
        transformations_str = "include file = " + self.root + "lib/geometry/transformations/" + transformations
        transformations_block = "{0}\n\t{1}\n{2}".format(start_transformations_delimiter,
                                                         transformations_str,
                                                         stop_transformations_delimiter)

        start_discovery_delimiter = ":start region discovery:"
        stop_discovery_delimiter = ":stop region discovery:"
        discovery_str = "action = " + discovery_action
        density_str = "density of random points (cm^-3) = {0}".format(density)
        boundary_str = "include file = {0}lib/geometry/sources/{1}/{2}/boundary.shape".format(self.root,
                                                                                              self.SOURCE_RADIONUCLIDE[
                                                                                                  self.source_model],
                                                                                              self.source_model)
        discovery_block = "{0}\n\t{1}\n\t{2}\n\t{3}\n{4}".format(start_discovery_delimiter,
                                                                 discovery_str,
                                                                 density_str,
                                                                 boundary_str,
                                                                 stop_discovery_delimiter)

        inscribed_geom_block = "{0}\n\t{1}\n\n\t{2}\n\n{3}".format(start_inscribed_geometry_delimiter,
                                                                   inscribed_geom_str,
                                                                   transformations_block,
                                                                   discovery_block,
                                                                   stop_inscribed_geometry_delimiter)

        autoenvelope_block = "{0}\n\t{1}\n\t{2}\n\t{3}\n\t{4}\n\n\t{5}\n{6}".format(start_geometry_delimiter,
                                                                                    name_str,
                                                                                    library_str,
                                                                                    type_str,
                                                                                    base_geometry_str,
                                                                                    inscribed_geom_block,
                                                                                    stop_geometry_delimiter)

        source_geom_str = "source geometries = seed"
        phantom_geom_str = "phantom geometries = phantom"
        switched_envelope_str = ""
        if superposition is True:
            switched_envelope_str = "source envelope geometry = phantom_and_seeds"

        simulation_geom_str = "simulation geometry = "



