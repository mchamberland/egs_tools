import sys
import os
import argparse
from os.path import join
import ct_tools.struct as cts
import ct_tools.ctdata as ctd
import ct_tools.artifact as cta
import ct_tools.resample as ctr
import ct_tools.ct_to_egsphant as cte
import brachy_dicom.reader as bdr


parser = argparse.ArgumentParser(description='Convert a DICOM CT dataset to the EGSnrc egsphant format.')

parser.add_argument('directory')
parser.add_argument('ctscheme')

parser.add_argument('--read_ctdata', dest='input_ctdata', nargs=1,
                    help='Read in a .ctdata file instead of a DICOM CT dataset. The .ctdata file is expected to be in '
                         'the directory provided as the first argument of the script.')

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

parser.add_argument('-v, --verbose', action='store_true',
                    help='Increase verbosity.')

args = parser.parse_args()


if not os.path.exists(args.directory):
    print('Directory {} does not exist.'.format(args.directory))
    sys.exit()

if not os.path.exists(args.ctscheme):
    print('CT conversion scheme file {} does not exist.'.format(args.ctscheme))
    sys.exit()


if args.directory == '.':
    current_directory = os.getcwd()
    base_name = os.path.basename(current_directory)
else:
    base_name = os.path.basename(args.directory)


if not args.input_ctdata:
    ctdata = ctd.get_ctdata_from_dicom(args.directory)
else:
    ctdata = ctd.read_ctdata(join(args.directory, args.input_ctdata))


if args.mar:
    search_radius, slices = 0, 0
    if len(args.mar) == 2:
        threshold, replacement = args.mar
    elif len(args.mar) == 4:
        threshold, replacement, search_radius, slices = args.mar
    else:
        raise Exception('Metal artifact reduction option expects the threshold and replacement values, and '
                        'optionally the search radius in mm and the number of z slices from the seed locations.')

    plan_filenames, plans = bdr.read_plan_files_in_directory(args.directory)
    if plans:
        seed_locations = bdr.get_seed_locations_in_mm(plans[0]) / 10
        artifact_reduction = cta.SimpleThresholdReplacement()
        if len(args.mar) == 2:
            artifact_reduction = cta.SimpleThresholdReplacement(threshold=threshold, replacement=replacement)
        else:
            artifact_reduction = cta.SimpleThresholdReplacement(threshold=threshold, replacement=replacement,
                                                                xy_search_in_mm=search_radius,
                                                                z_search_in_slices=slices)
        artifact_reduction.apply_str_to_seed_locations(ctdata, seed_locations)
    else:
        raise Exception('No DICOM RP plan files were found in directory')


if args.crop:
    ctdata = ctd.crop_ctdata_to_bounds(ctdata, args.crop)


if args.resample:
    nx, ny, nz, size_or_voxels = args.resample
    if not (size_or_voxels == 'size' or size_or_voxels == 'voxels'):
        raise Exception('Specify if the resampling is specified in cm (''size'') or in voxels (''voxels'').')
    if size_or_voxels == 'size':
        nx = float(nx)
        ny = float(ny)
        nz = float(nz)
    else:
        nx = int(nx)
        ny = int(ny)
        nz = int(nz)
    ctdata = ctr.resample_ctdata(ctdata, (nx, ny, nz), size_or_voxels)


contours = cts.get_contours_from_dicom(args.directory)


if args.write_ctdata:
    ctdata.write_to_file(join(args.directory, base_name))


ctconversion = cte.CTConversionToEGSphant(args.ctscheme, contours, is_verbose=args.verbose)
egsphant = ctconversion.convert_to_egsphant(ctdata)
egsphant.write_egsphant(join(args.directory, base_name))
