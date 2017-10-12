"""EGSinp

This module exports a single class that handles creation of egsinp files from templates for egs_brachy

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

    def __init__(self, filename="egsinp", mode="normal", egsphant=None, phsp_file=None):
        """Create a skeleton to build an egsinp file for egs_brachy"""
        self.root = os.path.expandvars("$EGS_HOME")
        self.run_mode = mode
        self.source_model = "OncoSeed_6711"
        self.filename = filename

        self.ncase = 1e6

        self.egsphant = egsphant
        self.box = None
        self.phantom = "10.0cmx10.0cmx10.0cm_2mm_xyz_water.geom"
        self.transformations = "single_seed_at_origin"

        self.volume_correction = "correct"

        self.phsp_file = phsp_file
        self.spectrum = "I125_NCRP_line"
        self.source_weights = ""

        self.score_tracklength_dose = "yes"
        self.score_edep = "no"
        self.muen_for_media = "WATER_0.998"

        self.reference_air_kerma_rate = 0.
        self.total_dwell_time = 0.
        self.dose_scaling_factor = 1.

        self.output_egsphant = "no"
        self.output_voxel_info = "no"

        self.times_to_reuse_recycled_particles = 1