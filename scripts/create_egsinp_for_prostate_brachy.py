#!/home/mchamber/anaconda3/bin/python
import sys
import os
from os.path import join
from shutil import copyfile
import argparse
import egsinp.egsinp as egsinp
from egsinp.egsinp import EGSinp
import dicom.reader as bdr
import ct_tools.ctdata as ctd


def read_seed_locations(filename):
    if not os.path.exists(filename):
        raise FileNotFoundError

    seed_locations_in_mm = []
    with open(filename) as file:
        lines = [line for line in file.readlines()]

        for line in lines:
            pos = [float(v)*10 for v in line.split()]
            seed_locations_in_mm.append(pos)

    return seed_locations_in_mm


parser = argparse.ArgumentParser(description='Create an egs_brachy input file based on a DICOM plan file and a '
                                             'CT-based egsphant.')

parser.add_argument('directory',
                    help='The directory where the DICOM RP file is located.')

parser.add_argument('-n', '--nhist',
                    help='The number of histories to run. This can (and probably should) be input in scientific '
                         'notation (e.g., 1e6).')

parser.add_argument('-c', '--copy_egsphant_to_egs_brachy_lib', dest='copy_egsphant', action='store_true',
                    help='Copy the egsphant to the egs_brachy/lib/geometry/phantoms/egsphant/ folder.')

parser.add_argument('-b', '--box_of_uniform_medium', dest='box', type=str,
                    help='A box of uniform medium (specified) in which the egsphant will be inscribed. The box matches'
                         'the dimensions and location of the CT dataset, by default. The medium needs to be specified.')

parser.add_argument('-s', '--seed_locations', dest='seed_locations', type=str,
                    help='A file containing the seed locations (x, y, z coordinates in cm; 1 per line).')

parser.add_argument('-r', '--rakr', dest='rakr', type=float,
                    help='The Reference Air-Kerma Rate (RAKR) of the seeds, in units of micro-gray per hour at 1 m.')

parser.add_argument('--volume-correction', dest='volcor', type=str,
                    help='The source volume correction option: "correct", "none", or "zero dose".')

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


if args.copy_egsphant:
    root_egsphant = join(os.path.expandvars("$EGS_HOME"), 'egs_brachy/lib/geometry/phantoms/egsphant')
    copyfile(join(args.directory, base_name + '.egsphant'), join(root_egsphant, base_name + '.egsphant'))

input_file = EGSinp(filename=join(args.directory, base_name), path_type='absolute')

plan_filenames, plans = bdr.read_plan_files_in_directory(args.directory)
if plans:
    plan = plans[0]
    if args.seed_locations:
        seed_locations = read_seed_locations(args.seed_locations)
    else:
        seed_locations = bdr.get_seed_locations_in_mm(plan)
    rakr = bdr.get_rakr_in_ugy_per_h_at_1m(plan)
    source_model = bdr.get_source_model(plan)

    egsinp.create_seed_transformations_file(seed_locations, filename=base_name, convert_to_cm=True,
                                            place_in_egs_brachy_lib=True)
    transformations = base_name + '.transf'
    media_list = 'PROSTATE_WW86, URETHRA_WW86, RECTUM_ICRP23, URINARY_BLADDER_EMPTY, P50C50, ' \
                 'CALCIFICATION_ICRU46, M_SOFT_TISSUE_ICRU46, CORTICAL_BONE_WW86'

    input_file.source_model = source_model
    if args.nhist:
        input_file.run_control(ncase=float(args.nhist))
    else:
        input_file.run_control()
    input_file.run_mode()
    input_file.media_definition()
    if args.box:
        ctdata = ctd.get_ctdata_from_dicom(args.directory)
        sx, sy, sz = ctdata.image_size_in_cm()
        cx, cy, cz = ctdata.image_centre_in_cm()
        # we pad the z-bounds slightly, to avoid overlapping boundaries with the egsphant
        egsinp.create_box_of_uniform_medium(base_name, (sx, sy, sz + 1.e-2), args.box, (cx, cy, cz))
        input_file.geometry(box=base_name, egsphant=base_name, transformations=transformations)
    else:
        input_file.geometry(egsphant=base_name, transformations=transformations)
    if args.volcor:
        input_file.source_volume_correction(correction=args.volcor)
    else:
        input_file.source_volume_correction(correction='zero dose')
    input_file.source(transformations=transformations)
    input_file.scoring_options(rakr=rakr, muen_file='brachy_xcom_1.5MeV_egsphant.muendat', muen_for_media=media_list)
    input_file.variance_reduction()
    input_file.transport_parameters()
    input_file.rng()
    input_file.close()

    root_egs_brachy = join(os.path.expandvars("$EGS_HOME"), 'egs_brachy/')
    copyfile(join(args.directory, base_name + '.egsinp'), join(root_egs_brachy, base_name + '.egsinp'))

else:
    raise Exception('No DICOM RT plan files were found in directory')


