#!/home/mchamber/anaconda3/bin/python
import sys
import os
import argparse
from os.path import join
import ct_tools.struct as cts
import ct_tools.ctdata as ctd
import ct_tools.artifact as cta
import ct_tools.resample as ctr
import ct_tools.ct_to_egsphant as cte
import dicom.reader as bdr
import dicom.rtdose as rtd


parser = argparse.ArgumentParser(description='Convert a DICOM CT dataset to the EGSnrc egsphant format.')

parser.add_argument('directory',
                    help='The directory where the DICOM (CT, RP, RS) or ctdata files are located.')
parser.add_argument('ctscheme',
                    help='The name of the ctscheme conversion file.')

parser.add_argument('-r', '--read_ctdata', dest='input_ctdata',
                    help='Read in a .ctdata file instead of a DICOM CT dataset. The .ctdata file is expected to be in '
                         'the directory provided as the first argument of the script.')

parser.add_argument('-a', '--apply_default_str', action='store_true',
                    help='Apply metal artifact reduction (Simple Threshold Reduction method, with default parameters.')

parser.add_argument('--apply_str', dest='mar', type=int, nargs='*',
                    help='Apply metal artifact reduction (Simple Threshold Reduction method). Optionally, specify the '
                         'threshold and replacement values and the search radius in mm and the number of z slices from '
                         'the seed locations.')

# TODO add option to specify box with 3 or 4 numbers (sizes in x, y, z, and optionally, centre of box)
parser.add_argument('-c', '--crop', type=float, nargs='*',
                    help='Crop the CT dataset to the (x1, x2, y1, y2, z1, z2) bounds specified (in cm). If a single '
                         'value n is passed, then the set is cropped to a cube of size n at the centre of the dataset')

parser.add_argument('--resample', nargs=4,
                    help='Resample the CT to the desired size (nx, ny, nz, type), where ''type'' specifies whether '
                         'the size is specified in cm (''size'') or in voxels (''voxels'').')

parser.add_argument('-m', '--match_dose_grid', dest='match', type=float, nargs=3,
                    help='The CT dataset will be resampled and cropped to match the dose grid within (tx, ty, tz).'
                         'For example, -m 1 -1.2 0 will resample the CT data to match the voxel size of the dose'
                         'grid. Then, it will crop the dataset to within 1 cm of the x-bounds of the dose grid;'
                         'it will crop an extra 1.2 cm along y (i.e., the cropped data will be shorter in y);'
                         'and it will crop to match exactly the z-bounds.')

parser.add_argument('-w', '--write_ctdata', action='store_true',
                    help='The CT data will be written to a text file before conversion to egsphant.')

parser.add_argument('-v', '--verbose', action='store_true',
                    help='Increase verbosity.')

args = parser.parse_args()


if not os.path.exists(args.directory):
    print('Directory {} does not exist.'.format(args.directory))
    sys.exit()

if not args.ctscheme.endswith('.ctscheme'):
    args.ctscheme += '.ctscheme'

if not os.path.exists(args.ctscheme):
    print('CT conversion scheme file {} does not exist.'.format(args.ctscheme))
    sys.exit()


if args.directory == '.':
    current_directory = os.getcwd()
    base_name = os.path.basename(current_directory)
else:
    base_name = os.path.basename(args.directory)

if args.verbose:
    print('Reading CT dataset...\n')
if not args.input_ctdata:
    ctdata = ctd.get_ctdata_from_dicom(args.directory)
else:
    ctdata = ctd.read_ctdata(join(args.directory, args.input_ctdata))


if args.apply_default_str:
    if args.verbose:
        print('Applying metal artifact reduction...\n')
    plan_filenames, plans = bdr.read_plan_files_in_directory(args.directory)
    if plans:
        seed_locations = bdr.get_seed_locations_in_mm(plans[0]) / 10
        artifact_reduction = cta.SimpleThresholdReplacement()
        artifact_reduction.apply_str_to_seed_locations(ctdata, seed_locations)
    else:
        raise Exception('No DICOM RP plan files were found in directory')


if args.mar:
    if args.verbose:
        print('Applying metal artifact reduction...')
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


if args.match:
    dose_filenames, doses = bdr.read_dose_files_in_directory(args.directory)
    if doses:
        dose = rtd.RTDoseInfo(doses[0])
        if args.verbose:
            print('Matching the dimensions and resolution to the dose grid...\n')
            print('Dose grid description:')
            dose.print_info()
        ctdata = ctr.resample_ctdata(ctdata, (dose.dx / 10, dose.dy / 10, dose.dz / 10), 'size')
        if args.verbose:
            print('Now cropping CT dataset...\n')
        dx, dy, dz = dose.grid_extents
        xi, xf = dx
        yi, yf = dy
        if dose.z_direction == -1:
            zf, zi = dz
        else:
            zi, zf = dz
        tx, ty, tz = args.match
        xi = xi / 10 - tx
        xf = xf / 10 + tx
        yi = yi / 10 - ty
        yf = yf / 10 + ty
        zi = zi / 10 - tz
        zf = zf / 10 + tz
        if any(args.match):
            ctdata = ctd.crop_ctdata_to_bounds(ctdata, (xi, xf, yi, yf, zi, zf))
        else:
            ctdata = ctd.crop_ctdata_to_bounds(ctdata, (xi, xf, yi, yf, zi, zf), include_partial_voxel=True)
        print('Modified CT dataset description:')
        ctdata.print_info()
        print()
    else:
        raise Exception('No DICOM RT dose files were found in directory')


if args.resample:
    if args.verbose:
        print('Resampling CT dataset...')
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


if args.crop:
    if args.verbose:
        print('Cropping CT dataset...')
    ctdata = ctd.crop_ctdata_to_bounds(ctdata, args.crop)


contours = cts.get_contours_from_dicom(args.directory)


if args.write_ctdata:
    if args.verbose:
        print('Writing CT data to file before conversion to egsphant starts...')
    ctdata.write_to_file(join(args.directory, base_name))


ctconversion = cte.CTConversionToEGSphant(args.ctscheme, contours, is_verbose=args.verbose)
egsphant = ctconversion.convert_to_egsphant(ctdata)
print('Writing egsphant...')
egsphant.write_egsphant(join(args.directory, base_name))
print('Done!')
