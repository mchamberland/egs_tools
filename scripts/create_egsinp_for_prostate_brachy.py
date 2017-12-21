#!/home/mchamber/venvs/py-3.6/bin/python
import sys
import os
import argparse
import egsinp.egsinp as egsinp
from egsinp.egsinp import EGSinp
import brachy_dicom.reader as bdr


parser = argparse.ArgumentParser(description='Create an egs_brachy input file based on a DICOM plan file and a '
                                             'CT-based egsphant.')

parser.add_argument('directory',
                    help='The directory where the DICOM RP file is located.')

parser.add_argument('-n', '--nhist',
                    help='The number of histories to run. This can (and probably should) be input in scientific '
                         'notation (e.g., 1e6).')

parser.add_argument('-v', '--verbose', action='store_true',
                    help='Increase verbosity.')


args = parser.parse_args()


if not os.path.exists(args.directory):
    print('Directory {} does not exist.'.format(args.directory))
    sys.exit()


if args.directory == '.':
    current_directory = os.getcwd()
    base_name = os.path.basename(current_directory)
else:
    base_name = os.path.basename(args.directory)


input_file = EGSinp(filename=base_name, path_type='absolute')

plan_filenames, plans = bdr.read_plan_files_in_directory(args.directory)
if plans:
    plan = plans[0]
    seed_locations = bdr.get_seed_locations_in_mm(plan)
    rakr = bdr.get_rakr_in_ugy_per_h_at_1m(plan)
    source_model = bdr.get_source_model(plan)

    egsinp.create_seed_transformations_file(seed_locations, filename=base_name, convert_to_cm=True,
                                            place_in_egs_brachy_lib=True)
    transformations = base_name + '.transf'
    media_list = 'PROSTATE_WW86, URETHRA_WW86, RECTUM_ICRP23, URINARY_BLADDER_EMPTY, URINARY_BLADDER_FULL, P50C50, ' \
                 'CALCIFICATION_ICRU46, M_SOFT_TISSUE_ICRU46, CORTICAL_BONE_WW86'

    input_file.source_model = source_model
    input_file.run_control(ncase=float(args.nhist))
    input_file.run_mode()
    input_file.media_definition()
    input_file.geometry(egsphant=base_name, transformations=transformations)
    input_file.source_volume_correction()
    input_file.source(transformations=transformations)
    input_file.scoring_options(rakr=rakr, muen_for_media=media_list)
    input_file.variance_reduction()
    input_file.transport_parameters()
    input_file.rng()

else:
    raise Exception('No DICOM RP plan files were found in directory')
