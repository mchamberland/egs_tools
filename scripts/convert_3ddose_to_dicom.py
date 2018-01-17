import sys
import os
import argparse
import dose_distribution.dose3d as dd

parser = argparse.ArgumentParser(description='Convert a 3ddose file to DICOM RT Dose format.')

parser.add_argument('dose',
                    help='The 3ddose file to convert.')

parser.add_argument('-d', '--dicom', type=str,
                    help='A DICOM RT dose file that will be used as a template. If no file is provided, a dummy '
                         'template DICOM RT dose file used. Note that if you try importing the dose in a TPS, it may be'
                         ' missing tags that are expected by TPS.')

parser.add_argument('-f', '--flipz', action='store_true',
                    help='3ddose files are stored in order of increasing z values. Often, DICOM CT and dose files are'
                         'stored in order of decreasing z values. This option will flip the order of the z slices of'
                         'of the 3ddose file when writing to the DICOM file.')

args = parser.parse_args()


dose = dd.DoseDistribution(args.dose)

template = None
if args.dicom:
    if os.path.exists(args.dicom):
        template = args.dicom
    else:
        print('DICOM file not found. Using a dummy template DICOM file instead.')

dicom = dd.write_3ddose_to_dicom(dose, template, args.flipz)
