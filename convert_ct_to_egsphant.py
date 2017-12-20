import sys
import os.path as path
import argparse
from os import listdir
from os.path import join, isfile, splitext
import ct_tools.struct as cts
import ct_tools.ctdata as ctd
import ct_tools.artifact as cta
import ct_tools.ct_to_egsphant as cte
import brachy_dicom.reader as bdr


parser = argparse.ArgumentParser(description='Convert a DICOM CT dataset to the EGSnrc egsphant format.')

parser.add_argument('directory')
parser.add_argument('ctscheme')

parser.add_argument('--apply_str', dest='mar', type=int, nargs='*',
                    help='Apply metal artifact reduction (Simple Threshold Reduction method). Optionally, specify the '
                         'threshold and replacement values and the search radius in mm and the number of z slices from '
                         'the seed locations.')

parser.add_argument('--crop', type=float, nargs='*',
                    help='Crop the CT dataset to the (x1, x2, y1, y2, z1, z2) bounds specified (in cm). If a single '
                         'value n is passed, then the set is cropped to a cube of size n at the centre of the dataset')

parser.add_argument('--resample', nargs=4,
                    help='Resample the CT to the desired size (nx, ny, nz, type), where ''type'' specifies whether '
                         'the size is specified in cm (''size'') or in voxels (''voxels''). Note that this is a very '
                         'computationally intensive and time consuming operation.')

parser.add_argument('--write_ctdata', action='store_true',
                    help='The CT data will be written to a text file before conversion to egsphant.')

parser.add_argument('--read_ctdata', dest='input_ctdata', nargs=1,
                    help='Read in a .ctdata file instead of a DICOM CT dataset.')

parser.add_argument('-v, --verbose', action='store_true',
                    help='Increase verbosity.')

args = parser.parse_args()

if not path.exists(args.directory):
    print('{} does not exist.'.format(args.directory))
    sys.exit()

if not args.ctscheme.endswith('.ctdata'):
    args.ctscheme += '.ctdata'

if not path.exists(args.ctscheme):
    print('{} does not exist.'.format(args.ctscheme))
    sys.exit()


ctdata = ctd.get_ctdata_from_dicom(args.directory)

contours = cts.get_contours_from_dicom(args.directory)

ctconversion = cte.CTConversionToEGSphant(args.ctscheme, contours, is_verbose=args.verbose)


files, plans = bdr.read_plan_files_in_directory('dicom')
seed_locations = bdr.get_seed_locations_in_mm(plans[0]) / 10
str = cta.SimpleThresholdReplacement()
str.apply_str_to_seed_locations(ctdata, seed_locations)
ctdata.write_to_file('1721177_MAR.ctdata')





ctconversion = cte.CTConversionToEGSphant('prostate_calc', contours, is_verbose=True)
egsphant = ctconversion.convert_to_egsphant(ctdata)
egsphant.write_egsphant('1721177.egsphant')