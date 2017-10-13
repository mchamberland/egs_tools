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

    t1 = "\n\t"
    t2 = "\n\t\t"
    t3 = "\n\t\t\t"
    t4 = "\n\t\t\t\t"

    def __init__(self, filename="egsinp", source_model="OncoSeed_6711", path_type="relative"):
        """Create a skeleton to build an egsinp file for egs_brachy"""
        self.root = ""
        self.source_model = source_model
        self.radionuclide = self.SOURCE_RADIONUCLIDE[self.source_model]
        self.airkerma_per_hist = self.SOURCE_AIR_KERMA_STRENGTH_PER_HISTORY[self.source_model]
        self.dwell_times_in_s = None
        self.total_dwell_time_in_s = 0.
        self.filename = filename
        if path_type == "absolute":
            self.root = os.path.expandvars("$EGS_HOME")
        self.input_file = open(filename + '.egsinp', 'a')

    def run_control(self, ncase=1e6, nbatch=1, nchunk=1, calculation='first', geometry_error_limit=2500,
                    egsdat_file_format='gzip'):
        start_delimiter = ":start run control:"
        stop_delimiter = ":stop run control:"
        ncase_str = "ncase = {:.0E}".format(ncase)
        nbatch_str = "nbatch = {}".format(nbatch)
        nchunk_str = "nchunk = {}".format(nchunk)
        calc_str = "calculation = {}".format(calculation)
        geom_error_str = "geometry error limit = {}".format(geometry_error_limit)
        egsdat_str = "egsdat file format = {}".format(egsdat_file_format)

        input_block = "\n{0}{t1}{1}{t1}{2}{t1}{3}{t1}{4}{t1}{5}{t1}{6}\n{7}\n".format(start_delimiter,
                                                                                      ncase_str,
                                                                                      nbatch_str,
                                                                                      nchunk_str,
                                                                                      calc_str,
                                                                                      geom_error_str,
                                                                                      egsdat_str,
                                                                                      stop_delimiter,
                                                                                      t1=self.t1)
        self.input_file.write(input_block)

    def run_mode(self, run_mode="normal", single_generator="yes"):
        start_delimiter = ":start run mode:"
        stop_delimiter = ":stop run mode:"
        mode_str = "run mode = {}".format(run_mode)
        gen_str = "single generator = {}".format(single_generator)
        input_block = "\n{0}{t1}{1}{t1}{2}\n{3}\n".format(start_delimiter,
                                                          mode_str,
                                                          gen_str,
                                                          stop_delimiter,
                                                          t1=self.t1)
        self.input_file.write(input_block)

    def transport_parameters(self, transport_file="low_energy_default"):
        input_block = "\ninclude file = {0}lib/transport/{1}\n".format(self.root,
                                                                       transport_file)
        self.input_file.write(input_block)

    def variance_reduction(self, do_recycling=True, recycling_number=1, deactivate_global_range_rejection=False,
                           global_range_rejection_max_total_energy=None):
        start_delimiter = ":start variance reduction:"
        stop_delimiter = ":stop variance reduction:"

        if do_recycling:
            recycling_str = "\t:start particle recycling:"
            recycling_str += "{1}times to reuse recycled particles = {0}".format(recycling_number, self.t2)
            recycling_str += "{}rotate recycled particles = yes".format(self.t2)
            recycling_str += "{}:stop particle recycling:".format(self.t1)
        else:
            recycling_str = ""

        range_rejection_str = ""
        if deactivate_global_range_rejection:
            range_rejection_str += "{}global range rejection = no".format(self.t1)
        if global_range_rejection_max_total_energy is not None:
            range_rejection_str += "{1}global range rejection max energy = {0}".format(
                global_range_rejection_max_total_energy, self.t1)

        input_block = "\n{0}\n{1}{2}\n{3}\n".format(start_delimiter,
                                                    recycling_str,
                                                    range_rejection_str,
                                                    stop_delimiter)
        self.input_file.write(input_block)

    def scoring_options(self, score_tracklength=True, score_edep=False, score_scatter=False,
                        muen_file="brachy_xcom_1.5MeV.muendat", muen_for_media="WATER_0.998", output_voxel_info=False,
                        rakr=None, dose_format="text"):
        start_delimiter = ":start scoring options:"
        stop_delimiter = ":stop scoring options:"

        if score_tracklength:
            track_str = "score tracklength dose = yes"
        else:
            track_str = "score tracklength dose = no"

        if score_edep:
            edep_str = "score energy deposition = yes"
        else:
            edep_str = "score energy deposition = no"

        if score_scatter:
            scatter_str = "score scatter dose = yes"
        else:
            scatter_str = "score scatter dose = no"

        muendat_str = "muen file = {0}lib/muen/{1}".format(self.root, muen_file)
        muen_media_str = "muen for media = {}".format(muen_for_media)

        dose_format_str = "dose file format = {}".format(dose_format)
        if output_voxel_info:
            voxel_str = "output voxel info files = yes"
            voxel_format_str = "voxel info file format = gzip"
        else:
            voxel_str = "output voxel info files = no"
            voxel_format_str = "voxel info file format = gzip"

        if rakr is not None:
            dose_scaling_factor = self.calculate_dose_scaling_factor(rakr)
            scaling_str = "dose scaling factor = ".format(dose_scaling_factor)
        else:
            scaling_str = "dose scaling factor = 1"

        input_block = "\n{0}{t1}{1}{t1}{2}{t1}{3}{t1}{4}{t1}{5}{t1}{6}{t1}{7}{t1}{8}{t1}{9}\n{10}\n".format(
            start_delimiter,
            track_str,
            edep_str,
            scatter_str,
            muendat_str,
            muen_media_str,
            voxel_str,
            voxel_format_str,
            scaling_str,
            dose_format_str,
            stop_delimiter,
            t1=self.t1)

        self.input_file.write(input_block)

    def media_definition(self, ae=0.512, ue=2.012, ap=0.001, up=1.500):
        start_delimiter = ":start media definition:"
        stop_delimiter = ":stop media definition:"
        ae_str = "AE = {}".format(ae)
        ue_str = "UE = {}".format(ue)
        ap_str = "AP = {}".format(ap)
        up_str = "UP = {}".format(up)
        material_str = "material data file = {0}lib/media/material.dat".format(self.root)
        input_block = "\n{0}{t1}{1}{t1}{2}{t1}{3}{t1}{4}{t1}{5}\n{6}\n".format(start_delimiter,
                                                                               ae_str,
                                                                               ue_str,
                                                                               ap_str,
                                                                               up_str,
                                                                               material_str,
                                                                               stop_delimiter,
                                                                               t1=self.t1)
        self.input_file.write(input_block)

    def source(self, source_type="isotropic", transformations="single_seed_at_origin", weights=None, phsp=None):
        start_delimiter = ":start source definition:\n\t:start source:"
        stop_delimiter = ":stop source definition:"
        stop_source_delimiter = ":stop source:"

        source_str = "name = {}".format(self.source_model)

        source_block = ""
        if source_type == "isotropic":
            library_str = "library = egs_isotropic_source"
            charge_str = "charge = 0"
            include_str = "include file = {0}lib/geometry/sources/{1}/{2}/{1}.shape".format(
                self.root,
                self.SOURCE_RADIONUCLIDE[self.source_model],
                self.source_model)

            start_spectrum_delimiter = ":start spectrum:"
            stop_spectrum_delimiter = ":stop spectrum:"
            spectrum_type_str = "type = tabulated spectrum"
            spectrum_file_str = "spectrum file = {0}lib/spectra/{1}.spectrum".format(
                self.root,
                self.RADIONUCLIDE_SPECTRUM[self.SOURCE_RADIONUCLIDE[self.source_model]])
            source_block = "{t2}{0}{t2}{1}{t2}{2}{t2}{3}{t2}{4}{t3}{5}{t3}{6}{t2}{7}{t1}{8}\n".format(
                source_str,
                library_str,
                charge_str,
                include_str,
                start_spectrum_delimiter,
                spectrum_type_str,
                spectrum_file_str,
                stop_spectrum_delimiter,
                stop_source_delimiter,
                t1=self.t1,
                t2=self.t2,
                t3=self.t3)
        elif phsp is not None:
            library_str = "library = eb_iaeaphsp_source"
            phsp_header_str = "header file = {0}lib/phsp/{1}.IAEAheader".format(self.root, phsp)
            source_block = "{t2}{0}{t2}{1}{t2}{2}{t1}{3}\n".format(source_str,
                                                                   library_str,
                                                                   phsp_header_str,
                                                                   stop_source_delimiter,
                                                                   t2=self.t2,
                                                                   t1=self.t1)

        start_transformations_delimiter = ":start transformations:"
        stop_transformations_delimiter = ":stop transformations:"
        transformations_str = "include file = {0}lib/geometry/transformations/{1}".format(self.root, transformations)
        transformations_block = "{t1}{0}{t2}{1}{t1}{2}\n".format(start_transformations_delimiter,
                                                                 transformations_str,
                                                                 stop_transformations_delimiter,
                                                                 t1=self.t1,
                                                                 t2=self.t2)

        weights_str = ""
        if weights is not None:
            the_weights = " ".join(str(w) for w in weights)
            weights_str = "source weights = {}".format(the_weights)
        elif self.dwell_times_in_s is not None:
            weights = self.get_source_weights_from_dwell_times()
            the_weights = " ".join(str(w) for w in weights)
            weights_str = "source weights = {}".format(the_weights)

        simulation_source_str = "simulation source = {}".format(self.source_model)

        input_block = "\n{0}{1}{2}{t1}{3}{t1}{4}\n{5}\n".format(start_delimiter,
                                                                source_block,
                                                                transformations_block,
                                                                weights_str,
                                                                simulation_source_str,
                                                                stop_delimiter,
                                                                t1=self.t1)

        self.input_file.write(input_block)

    def source_volume_correction(self, correction="correct", density=1e8):
        start_delimiter = ":start volume correction:\n\t:start source volume correction:"
        stop_delimiter = ":stop source volume correction:\n:stop volume correction:"

        correction_str = "correction type = {0}".format(correction)
        density_str = "density of random points (cm^-3) = {0}".format(density)
        boundary_str = "include file = {0}lib/geometry/sources/{1}/{2}/boundary.shape"\
            .format(self.root,
                    self.SOURCE_RADIONUCLIDE[self.source_model],
                    self.source_model)

        input_block = "\n{0}{t2}{1}{t2}{2}{t2}{3}{t1}{4}\n".format(start_delimiter,
                                                                   correction_str,
                                                                   density_str,
                                                                   boundary_str,
                                                                   stop_delimiter,
                                                                   t1=self.t1,
                                                                   t2=self.t2)
        self.input_file.write(input_block)

    def extra_volume_correction(self, correction="correct", density=1e8, boundary="boundary"):
        start_delimiter = ":start volume correction:\n\t:start extra volume correction:"
        stop_delimiter = ":stop extra volume correction:\n:stop volume correction:"

        correction_str = "correction type = {0}".format(correction)
        density_str = "density of random points (cm^-3) = {0}".format(density)
        boundary_str = "include file = {}.shape".format(boundary)

        input_block = "\n{0}{t2}{1}{t2}{2}{t2}{3}{t1}{4}\n".format(start_delimiter,
                                                                   correction_str,
                                                                   density_str,
                                                                   boundary_str,
                                                                   stop_delimiter,
                                                                   t1=self.t1,
                                                                   t2=self.t2)
        self.input_file.write(input_block)

    def volume_correction_from_file(self, phantoms=None, volcor=None):
        start_delimiter = ":start volume correction:\n\t:start volume correction from file:"
        stop_delimiter = ":stop volume correction from file:\n:stop volume correction:"
        phantom_str = ""
        for p, v in zip(phantoms, volcor):
            phantom_str += "{t2}phantom file = {0} {1}".format(p, v, t2=self.t2)

        input_block = "\n{0}{1}{t1}{2}\n".format(start_delimiter,
                                                 phantom_str,
                                                 stop_delimiter,
                                                 t1=self.t1)
        self.input_file.write(input_block)

    def geometry(self, box=None, phantom=None, egsphant=None, transformations="single_seed_at_origin",
                 discovery_action="discover", density=1e8, superposition=False):
        start_delimiter = ":start geometry definition:"
        stop_delimiter = ":stop geometry definition:"
        start_geometry_delimiter = ":start geometry:"
        stop_geometry_delimiter = ":stop geometry:"

        box_block = ""
        if box is not None:
            name_str = "name = box"
            library_str = "library = egs_glib"
            include_str = "include file = {0}lib/geometry/phantoms/{1}.geom".format(self.root, box)
            box_block = "{t1}{0}{t2}{1}{t2}{2}{t2}{3}{t1}{4}\n".format(start_geometry_delimiter,
                                                                       name_str,
                                                                       library_str,
                                                                       include_str,
                                                                       stop_geometry_delimiter,
                                                                       t1=self.t1,
                                                                       t2=self.t2)
        phantom_block = ""
        if phantom is not None:
            name_str = "name = phantom"
            library_str = "library = egs_glib"
            include_str = "include file = {0}lib/geometry/phantoms/{1}.geom".format(self.root, phantom)
            phantom_block = "{t1}{0}{t2}{1}{t2}{2}{t2}{3}{t1}{4}\n".format(start_geometry_delimiter,
                                                                           name_str,
                                                                           library_str,
                                                                           include_str,
                                                                           stop_geometry_delimiter,
                                                                           t1=self.t1,
                                                                           t2=self.t2)
        egsphant_block = ""
        if egsphant is not None:
            name_str = "name = phantom"
            library_str = "library = egs_glib"
            type_str = "type = egsphant"
            egsphant_str = "egsphant file = {0}lib/geometry/phantoms/egsphant/{1}.egsphant".format(self.root, egsphant)
            egsphant_block = "{t1}{0}{t2}{1}{t2}{2}{t2}{3}{t1}{4}\n".format(start_geometry_delimiter,
                                                                            name_str,
                                                                            library_str,
                                                                            type_str,
                                                                            egsphant_str,
                                                                            stop_geometry_delimiter,
                                                                            t1=self.t1,
                                                                            t2=self.t2)

        name_str = "name = seed"
        library_str = "library = egs_glib"
        include_str = "include file = {0}lib/geometry/sources/{1}/{2}/{2}.shape".format(
            self.root, self.SOURCE_RADIONUCLIDE[self.source_model], self.source_model)
        seed_block = "{t1}{0}{t2}{1}{t2}{2}{t2}{3}{t1}{4}\n".format(start_geometry_delimiter,
                                                                    name_str,
                                                                    library_str,
                                                                    include_str,
                                                                    stop_geometry_delimiter,
                                                                    t1=self.t1,
                                                                    t2=self.t2)

        name_str = "name = phantom_and_seeds"
        library_str = "library = egs_autoenvelope"
        type_str = "type = EGS_AEnvelope"
        if superposition is True:
            type_str = "type = EGS_ASwitchedEnvelope"
        base_geometry_str = "base geometry = phantom"

        start_inscribed_geometry_delimiter = ":start inscribed geometry:"
        stop_inscribed_geometry_delimiter = ":stop inscribed geometry:"
        inscribed_geom_str = "inscribed geometry name = seed"

        start_transformations_delimiter = ":start transformations:"
        stop_transformations_delimiter = ":stop transformations:"
        transformations_str = "include file = {0}lib/geometry/transformations/{1}".format(self.root, transformations)
        transformations_block = "{t3}{0}{t4}{1}{t3}{2}\n".format(start_transformations_delimiter,
                                                                 transformations_str,
                                                                 stop_transformations_delimiter,
                                                                 t3=self.t3,
                                                                 t4=self.t4)

        start_discovery_delimiter = ":start region discovery:"
        stop_discovery_delimiter = ":stop region discovery:"
        discovery_str = "action = {}".format(discovery_action)
        density_str = "density of random points (cm^-3) = {:.0E}".format(density)
        boundary_str = "include file = {0}lib/geometry/sources/{1}/{2}/boundary.shape".format(self.root,
                                                                                              self.SOURCE_RADIONUCLIDE[
                                                                                                  self.source_model],
                                                                                              self.source_model)
        discovery_block = "{t3}{0}{t4}{1}{t4}{2}{t4}{3}{t3}{4}".format(start_discovery_delimiter,
                                                                       discovery_str,
                                                                       density_str,
                                                                       boundary_str,
                                                                       stop_discovery_delimiter,
                                                                       t3=self.t3,
                                                                       t4=self.t4)

        inscribed_geom_block = "{t2}{0}{t3}{1}{2}{3}\n".format(start_inscribed_geometry_delimiter,
                                                               inscribed_geom_str,
                                                               transformations_block,
                                                               discovery_block,
                                                               stop_inscribed_geometry_delimiter,
                                                               t2=self.t2,
                                                               t3=self.t3)

        autoenvelope_block = "{t1}{0}{t2}{1}{t2}{2}{t2}{3}{t2}{4}{5}{t1}{6}\n".format(start_geometry_delimiter,
                                                                                      name_str,
                                                                                      library_str,
                                                                                      type_str,
                                                                                      base_geometry_str,
                                                                                      inscribed_geom_block,
                                                                                      stop_geometry_delimiter,
                                                                                      t1=self.t1,
                                                                                      t2=self.t2)
        envelope_block = ""
        if box is not None:
            name_str = "name = final"
            library_str = "library = egs_genvelope"
            base_geometry_str = "base geometry = box"
            inscribed_geom_str = "inscribed geometries = phantom_and_seeds"
            envelope_block = "{t1}{0}{t2}{1}{t2}{2}{t2}{3}{t2}{4}{t1}{5}\n".format(start_geometry_delimiter,
                                                                                   name_str,
                                                                                   library_str,
                                                                                   base_geometry_str,
                                                                                   inscribed_geom_str,
                                                                                   stop_geometry_delimiter,
                                                                                   t1=self.t1,
                                                                                   t2=self.t2)
            simulation_geom_str = "simulation geometry = final"
        else:
            simulation_geom_str = "simulation geometry = phantom_and_seeds"

        source_geom_str = "source geometries = seed"
        phantom_geom_str = "phantom geometries = phantom"
        switched_envelope_str = ""
        if superposition is True:
            switched_envelope_str = "source envelope geometry = phantom_and_seeds"

        input_block = "\n{0}{1}{2}{3}{4}{5}{6}{t1}{7}{t1}{8}{t1}{9}{t1}{10}\n{11}\n".format(start_delimiter,
                                                                                            box_block,
                                                                                            phantom_block,
                                                                                            egsphant_block,
                                                                                            seed_block,
                                                                                            autoenvelope_block,
                                                                                            envelope_block,
                                                                                            source_geom_str,
                                                                                            phantom_geom_str,
                                                                                            switched_envelope_str,
                                                                                            simulation_geom_str,
                                                                                            stop_delimiter,
                                                                                            t1=self.t1)
        self.input_file.write(input_block)

    def rng(self, seed1=29, seed2=33):
        start_delimiter = ":start rng definition:"
        stop_delimiter = ":stop rng definition:"
        seeds_str = "initial seeds = {0} {1}".format(seed1, seed2)
        input_block = "\n{0}{t1}{1}\n{2}\n".format(start_delimiter,
                                                   seeds_str,
                                                   stop_delimiter,
                                                   t1=self.t1)
        self.input_file.write(input_block)

    def calculate_dose_scaling_factor(self, rakr=1.):
        if self.SOURCE_RADIONUCLIDE[self.source_model].endswith("LDR"):
            lifetime = self.MEAN_LIFETIME_IN_HOURS[self.radionuclide]
            dose_scaling_factor = rakr * lifetime / self.airkerma_per_hist / 100

        elif self.SOURCE_RADIONUCLIDE[self.source_model].endswith("HDR"):
            dose_scaling_factor = rakr * self.total_dwell_time_in_s * max(self.get_source_weights_from_dwell_times())\
                                  / self.airkerma_per_hist / 3600. / 100.

        else:
            dose_scaling_factor = 1.

        return dose_scaling_factor

    def get_source_weights_from_dwell_times(self):
        weights = [t / self.dwell_times_in_s for t in self.dwell_times_in_s]
        return weights

    def set_dwell_times_in_s(self, dwell_times=None):
        self.dwell_times_in_s = dwell_times
        for t in self.dwell_times_in_s:
            self.total_dwell_time_in_s += t

    def close(self):
        self.input_file.close()

    def open(self):
        self.input_file = open(self.filename + '.egsinp', 'a')
