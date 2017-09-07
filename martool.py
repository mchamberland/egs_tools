#!/usr/bin/env python
import sys
import pydicom
import numpy
import os
from os import walk
from scipy import spatial
import matplotlib
import operator
from matplotlib.path import Path
import matplotlib.pyplot as plt
matplotlib.use('TkAgg')


def find_nearest(array, value):
    idx = (numpy.abs(array - value)).argmin()
    return array[idx]


#########
# File Structure:
# /...directory.../Patient[#]/CT/CT*.dicom
# /...directory.../Patient[#]/(RTP,RTS,RTD).dicom
#########

labelID = 'PC'

PatientNumbers = 'all'

if PatientNumbers == 'all':
    pts = os.listdir('/Users/Marc/clrp/prostate_cases_Mehan/patients')
    PatientNumbers = []

    for i in pts:
        if i != '.DS_Store':
            PatientNumbers.append(i.split('Patient')[1])
else:
    PatientNumbers = PatientNumbers.split(',')

for PatientNumber in PatientNumbers:

    print('***Patient ' + PatientNumber + '***')

    FilePath = '/Users/Marc/clrp/prostate_cases_Mehan/patients/Patient' + PatientNumber + '/'
    CTFilePath = FilePath + 'CT/'

    #
    # Open CT DICOM Files#
    #

    ctfilenamelist = []
    for (dirpath, dirnames, filenames) in walk(CTFilePath):
        ctfilenamelist.extend(filenames)
        break

    ctslice = []
    header = ctfilenamelist[1].split('.')
    del header[-1]
    del header[-1]
    headername = ".".join(header)
    for i in range(0, len(ctfilenamelist) - 1):  # reversed(range(0,len(ctfilenamelist))):
        if i != '.DS_Store':
            ctslice.append(pydicom.read_file(CTFilePath + headername + "." + str(i) + '.dcm'))  # Open them from xx.dcm to 0.dcm

    #######################
    # Open Dose, Struct, Plan DICOM Files#
    #######################

    filenamelist = []
    for (dirpath, dirnames, filenames) in walk(FilePath):
        filenamelist.extend(filenames)
        break

    for filename in filenamelist:
        if filename[0:3] == 'RTD' or filename[0:2] == 'RD':
            rd = pydicom.read_file(FilePath + filename)
        if filename[0:3] == 'RTS' or filename[0:2] == 'RS':
            rs = pydicom.read_file(FilePath + filename)
        if filename[0:3] == 'RTP' or filename[0:2] == 'RP':
            rp = pydicom.read_file(FilePath + filename)

    # Important Info#
    m = [ctslice[0].Columns, ctslice[0].Rows, len(ctslice)]

    n = [ctslice[0].PixelSpacing[0] / 10, ctslice[0].PixelSpacing[1] / 10,
         abs((ctslice[1].ImagePositionPatient[2] - ctslice[0].ImagePositionPatient[2]) / 10)]

    ############
    # STR MAR Method#
    ############

    # Set some values here
    tolerance = 0.025
    maxdistancexy = n[0] * 10
    maxdistancez = n[2] * 2

    # Some options
    STR = 243  # Threshold
    REPL = 25  # Replacement value

    xyzbounds = numpy.zeros((m[2], m[1]), dtype=object)
    xstart = ctslice[0].ImagePositionPatient[0] / 10
    ystart = ctslice[0].ImagePositionPatient[1] / 10
    zstart = ctslice[0].ImagePositionPatient[2] / 10

    xbounds = [ctslice[0].ImagePositionPatient[0] / 10]
    for i in range(0, m[0]):
        xbounds.append((xbounds[i] + n[0]))

    ybounds = [ctslice[0].ImagePositionPatient[1] / 10]
    for i in range(0, m[1]):
        ybounds.append((ybounds[i] + n[1]))

    zbounds = [ctslice[0].ImagePositionPatient[2] / 10]
    for i in range(0, m[2]):
        zbounds.append(round((zbounds[i] - n[2]), 1))

    with open('bounds.txt', 'a') as f:
        counter = 0
        for x in xbounds:
            if counter % 3 == 0:
                f.write('\n')
            f.write(str(x) + '\t')
            counter += 1

        counter = 0
        for y in ybounds:
            if counter % 3 == 0:
                f.write('\n')
            f.write(str(y) + '\t')
            counter += 1

        counter = 0
        for z in reversed(zbounds):
            if counter % 3 == 0:
                f.write('\n')
            f.write(str(z))
            counter += 1

    seedlocation = rp.ApplicationSetupSequence
    xnearest = [[] for x in range(0, len(seedlocation))]
    ynearest = [[] for y in range(0, len(seedlocation))]
    znearest = [[] for z in range(0, len(seedlocation))]
    locations = [[] for z in range(0, m[2])]

    for i in range(0, len(seedlocation)):
        xes = list(xbounds)
        yes = list(ybounds)
        zes = list(zbounds)

        x = round(seedlocation[i].ChannelSequence[0].BrachyControlPointSequence[0].ControlPoint3DPosition[0] / 10, 4)
        y = round(seedlocation[i].ChannelSequence[0].BrachyControlPointSequence[0].ControlPoint3DPosition[1] / 10, 4)
        z = round(seedlocation[i].ChannelSequence[0].BrachyControlPointSequence[0].ControlPoint3DPosition[2] / 10, 4)
        locations[abs(int(round(z / n[2], 0)))].append((x, y))

        dis = 0
        while dis < maxdistancexy:
            anchor = x  # seed_x
            xboundsnp = numpy.asarray(xes)
            point = find_nearest(xboundsnp, anchor)
            # xnearest.append(point)
            xnearest[i].append(xes.index(point))
            xes[xes.index(point)] = 999
            # xbounds.remove(point)
            dis = point - anchor

        dis = 0
        while dis < maxdistancexy:
            anchor = y  # seed_y
            yboundsnp = numpy.asarray(yes)
            point = find_nearest(yboundsnp, anchor)
            # ynearest.append(point)
            ynearest[i].append(yes.index(point))
            yes[yes.index(point)] = 999
            # yes.remove(point)
            dis = point - anchor

        dis = 0
        while dis < maxdistancez:
            anchor = z  # seed_y
            zboundsnp = numpy.asarray(zes)
            point = find_nearest(zboundsnp, anchor)
            # znearest.append(point)
            znearest[i].append(zes.index(point))
            zes[zes.index(point)] = 999
            # zes.remove(point)
            dis = point - anchor

    counter = 0
    print("Starting STR for Threshold = " + str(STR) + " and Replacement = " + str(REPL))
    for i in range(0, len(xnearest)):
        counter += 1
        sys.stdout.write("\r" + str(round(counter / len(xnearest) * 100, 2)) + "%")
        sys.stdout.flush()
        for z in znearest[i]:
            for y in ynearest[i]:
                for x in xnearest[i]:
                    try:
                        if ctslice[z].pixel_array[y][x] > STR:
                            ctslice[z].pixel_array[y][x] = REPL
                    except:
                        pass
    print("...done!")

    with open('test.CTdata', 'a') as f:
        for slice in reversed(range(0, len(ctslice))):
            numpy.savetxt(f, ctslice[slice].pixel_array, fmt='%i', footer='', comments='')

    print("All done!")


    ##############
    # Contour Interpolation#
    ##############

    # Checking contour bounds

    contourslices = [[[] for xx in range(int(round(n[2] * m[2] / 0.2, 0)))] for x in
                     range(0, len(rs.ROIContourSequence))]
    Label = []
    for ROI in range(0, len(rs.ROIContourSequence)):
        Label.append(rs.RTROIObservationsSequence[ROI].ROIObservationLabel)
        for Sequence in range(0, len(rs.ROIContourSequence[ROI].ContourSequence)):
            slicenumber2 = abs((rs.ROIContourSequence[ROI].ContourSequence[Sequence].ContourData[2] / 10 / n[2]))
            slicenumber = int(
                round(abs((rs.ROIContourSequence[ROI].ContourSequence[Sequence].ContourData[2] / 10 / n[2])), 0))
            for i in range(0, len(rs.ROIContourSequence[ROI].ContourSequence[Sequence].ContourData)):
                if i % 3 == 0:
                    contourslices[ROI][slicenumber].append((round(
                        rs.ROIContourSequence[ROI].ContourSequence[Sequence].ContourData[i] / 10, 5), round(
                        rs.ROIContourSequence[ROI].ContourSequence[Sequence].ContourData[i + 1] / 10, 5)))

    # Interpolate for slices with no contours
    print("Interpolating contours...")
    for ROI in range(0, len(rs.ROIContourSequence)):
        sys.stdout.write("..." + Label[ROI])
        sys.stdout.flush()
        exist = []
        for i in range(len(contourslices[ROI])):
            if contourslices[ROI][i]:
                exist.append(i)
        existpoints = [i for i in range(exist[0], exist[-1] + 1)]
        for i in existpoints:
            if not contourslices[ROI][i]:
                j = i - 1
                k = i + 1
                while not contourslices[ROI][j]:
                    j -= 1
                while not contourslices[ROI][k]:
                    k += 1
                dist, idx = spatial.cKDTree(contourslices[ROI][j]).query(contourslices[ROI][k])
                kcount = 0
                for element in idx:
                    jpoint = contourslices[ROI][j][element]
                    kpoint = contourslices[ROI][k][kcount]
                    weight = 1 / (k - j)
                    jpoint_weighted = tuple([(1 - weight) * x for x in jpoint])
                    kpoint_weighted = tuple([weight * x for x in kpoint])
                    inter = tuple(map(operator.add, jpoint_weighted, kpoint_weighted))
                    contourslices[ROI][i].append((round(inter[0], 5), round(inter[1], 5)))
                    kcount += 1
    print("...done!")

    ##########################
    # Contrast Bladder Checker#
    ##########################

    # Create Contour Mask# Needed for egsphant, so don't remove!
    print("Creating contour mask...")
    xybounds = []
    for y in ybounds:
        for x in xbounds:
            xybounds.append((x, y))
    cmap = [[[] for z in range(0, len(ybounds))] for ROI in range(0, 4)]  # range(0,len(contourslices))]
    for i in range(0, len(cmap)):
        sys.stdout.write("\r" + str(round((i + 1) / len(cmap) * 100, 2)) + "%")
        sys.stdout.flush()
        for j in range(0, m[2]):
            try:
                cmap[i][j] = Path(contourslices[i][j]).contains_points(xybounds)
            except:
                cmap[i][j] = [False for x in range(0, len(xbounds) * len(ybounds))]

    print('...done!')

    roibladder = -1
    for ROI in Label:
        if ROI.lower() == 'bladder':
            roibladder = Label.index(ROI)

    roiprostate = -1
    for ROI in Label:
        if ROI.lower() == 'target' or ROI.lower() == 'prostate':
            roiprostate = Label.index(ROI)

    roirectum = -1
    for ROI in Label:
        if ROI.lower() == 'rectum':
            roirectum = Label.index(ROI)

    roiurethra = -1
    for ROI in Label:
        if ROI.lower() == 'urethra':
            roiurethra = Label.index(ROI)

    if roibladder != -1:
        bladderctnumber = []
        for bladdercheck in contourslices[roibladder]:
            if bladdercheck:
                sl = contourslices[roibladder].index(bladdercheck)
                for y in range(0, len(ctslice[sl].pixel_array)):
                    for x in range(0, len(ctslice[sl].pixel_array[y])):
                        if cmap[roibladder][sl][x + y * len(ybounds)]:
                            bladderctnumber.append(ctslice[sl].pixel_array[y][x])
                break

        bladderavg = sum(bladderctnumber) / len(bladderctnumber)
        bladderflag = 0
        if bladderavg > 100:  # if the CT HU is greater than a threshold value, we know this patient has contrast in
            # their bladder
            bladderflag = 1
            print('**Warning**\tPatient' + PatientNumber + ' has contrast present in their bladder')

    if bladderflag == 1:
        print("Creating feathered contour mask...")
        xybounds = []
        for y in ybounds:
            for x in xbounds:
                xybounds.append((x, y))
        cmap = [[[] for z in range(0, len(ybounds))] for ROI in range(0, len(contourslices))]
        for i in range(0, len(cmap)):
            sys.stdout.write("\r" + str(round((i + 1) / len(cmap) * 100, 2)) + "%")
            sys.stdout.flush()
            for j in range(0, m[2]):
                try:
                    cmap[i][j] = Path(contourslices[i][j]).contains_points(xybounds, radius=-1.0)
                except:
                    cmap[i][j] = [False for x in range(0, len(xbounds) * len(ybounds))]
        print('...done!')

        print("Replacing contrast bladder...")
        for z in range(sl, len(ctslice)):
            sys.stdout.write("\r" + str(round((z - sl + 1) / (len(ctslice) - sl) * 100, 2)) + "%")
            sys.stdout.flush()
            for y in range(0, len(ctslice[z].pixel_array)):
                for x in range(0, len(ctslice[z].pixel_array[y])):
                    ex = ctslice[z].pixel_array[y][x]
                    if cmap[roibladder][z][x + y * len(ybounds)] == True and ex > 50:
                        ctslice[z].pixel_array[y][x] = 8
        print('...done!')

        print("Creating normal contour mask...")
        xybounds = []
        for y in ybounds:
            for x in xbounds:
                xybounds.append((x, y))
        cmap = [[[] for z in range(0, len(ybounds))] for ROI in range(0, 4)]  # range(0,len(contourslices))]
        for i in range(0, len(cmap)):
            sys.stdout.write("\r" + str(round((i + 1) / len(cmap) * 100, 2)) + "%")
            sys.stdout.flush()
            for j in range(0, m[2]):
                try:
                    cmap[i][j] = Path(contourslices[i][j]).contains_points(xybounds)
                except:
                    cmap[i][j] = [False for x in range(0, len(xbounds) * len(ybounds))]

        print('...done!')

    #
    # Write CT images#
    #

    colors = ['.g', '.y', '.b', '.r', '.w', '.c']
    if not os.path.exists(FilePath + 'MAR_img_' + labelID + '/'):
        os.makedirs(FilePath + 'MAR_img_' + labelID + '/')
    for i in range(0, len(ctslice)):
        sys.stdout.write('\r' + "...slice " + str(i + 1) + "/" + str(len(ctslice)))
        sys.stdout.flush()
        plt.figure()
        fig = plt.imshow(ctslice[i].pixel_array, cmap=plt.cm.gray)
        for ROI in range(0, len(contourslices)):
            if contourslices[ROI][i]:
                xpoints = [(x[0] / n[0] + (m[0] / 2)) for x in contourslices[ROI][i]]
                ypoints = [(y[1] / n[1] + (m[1] / 2)) for y in contourslices[ROI][i]]
                plt.plot(xpoints, ypoints, colors[ROI])
        if locations[i]:
            xpoints = [(x[0] / n[0] + (m[0] / 2)) for x in locations[i]]
            ypoints = [(y[1] / n[1] + (m[1] / 2)) for y in locations[i]]
            plt.plot(xpoints, ypoints, '*k')
        plt.xlim(m[0], 0)
        plt.ylim(m[1], 0)
        plt.axis('off')
        fig.axes.get_xaxis().set_visible(False)
        fig.axes.get_yaxis().set_visible(False)
        string = 'Patient' + PatientNumber + ' :: Slice ' + str(i + 1) + '/' + str(len(ctslice))
        plt.text(m[0] - 10, 20, string, color='white')
        plt.tight_layout()
        plt.savefig(FilePath + 'MAR_img_' + labelID + '/' + 'Slice_' + "%03d" % (i,) + '.png', bbox_inches='tight',
                    pad_inches=0)
        plt.close()
    print('...done!')

    #
    # Make EGSPHANT
    #

    Media = ['PROSTATE_WW86', 'CALCIFICATION_ICRU46', 'URINARY_BLADDER_EMPTY', 'AIR_TG43', 'RECTUM_ICRP23',
             'M_SOFT_TISSUE_ICRU46', 'CORTICAL_BONE_WW86',
             'P50C50']  # Need to create multiple tissue assignment schemes
    MediaDict = {'Target': 1, 'Rectum': 5, 'Bladder': 3, 'Urethra': 1}
    BigBore = [(-1024, 0.010), (-706, 0.280), (-535, 0.400), (-99, 0.900), (10, 1.000), (208, 1.090), (467, 1.280),
               (1234, 1.690), (3520, 2.600), (5750, 3.500), (8000, 4.380), (16000, 7.500)]  # Electron Densities

    BigBore_x = []
    BigBore_y = []
    for i in BigBore:
        BigBore_x.append(i[0])
        BigBore_y.append(-0.1746 + (1.176 * i[1]))  # Conversion to mass densities

    phantpath = FilePath + 'Patient' + PatientNumber + '_' + labelID + '.egsphant'
    egsphant = open(phantpath, 'w')

    # Header Data#
    egsphant.write(str(len(Media)) + '\n')  # Number of Media
    for tissue in Media:  # List of Tissues
        egsphant.write(tissue + '\n')

    egsphant.write(' ' * 2 + '0.25')
    for i in range(0, len(Media) - 1):
        egsphant.write(' ' * 7 + '0.25')
    egsphant.write('\n')
    for i in m:
        if i < 100:
            egsphant.write(' ')
        egsphant.write(' ' * 2 + str(i))

    counter = 0
    for x in xbounds:
        if counter % 3 == 0:
            egsphant.write('\n')
            egsphant.write(' ' * 3)
        egsphant.write("{:.5f}".format(x))
        if counter % 3 == 0 or counter % 3 == 1:
            egsphant.write(' ' * 10)
        if counter % 3 == 2:
            egsphant.write(' ' * 7)
        counter += 1

    counter = 0
    for y in ybounds:
        if counter % 3 == 0:
            egsphant.write('\n')
            egsphant.write(' ' * 3)
        egsphant.write("{:.5f}".format(y))
        if counter % 3 == 0 or counter % 3 == 1:
            egsphant.write(' ' * 10)
        if counter % 3 == 2:
            egsphant.write(' ' * 7)
        counter += 1

    counter = 0
    for z in reversed(zbounds):
        if counter % 3 == 0:
            egsphant.write('\n')
            egsphant.write(' ' * 3)
        egsphant.write("{:.5f}".format(z))
        if counter % 3 == 0 or counter % 3 == 1:
            egsphant.write(' ' * 10)
        if counter % 3 == 2:
            egsphant.write(' ' * 7)
        counter += 1
    egsphant.write('\n')

    print('Writing voxel tissue assignments...')
    for z in reversed(range(0, len(ctslice))):
        sys.stdout.write("\r" + str(round(abs(z - m[2]) / m[2] * 100, 2)) + "%")
        sys.stdout.flush()
        for y in range(0, len(ctslice[z].pixel_array)):
            for x in reversed(range(0, len(ctslice[z].pixel_array[y]))):
                ex = ctslice[z].pixel_array[y][x]
                if cmap[roiurethra][z][x + y * len(xbounds)]:
                    egsphant.write(str(MediaDict['Urethra']))
                    continue
                elif cmap[roiprostate][z][x + y * len(xbounds)]:
                    if ex < 244:
                        egsphant.write(str(MediaDict['Target']))
                        continue
                    if 243 < ex < 399:
                        egsphant.write(str(Media.index('P50C50') + 1))
                        continue
                    if ex > 398:
                        egsphant.write(str(Media.index('CALCIFICATION_ICRU46') + 1))
                        continue
                elif cmap[roirectum][z][x + y * len(xbounds)]:
                    egsphant.write(str(MediaDict['Rectum']))
                    continue
                elif cmap[roibladder][z][x + y * len(xbounds)]:
                    egsphant.write(str(MediaDict['Bladder']))
                    continue
                elif ex < 244:
                    egsphant.write(str(Media.index('M_SOFT_TISSUE_ICRU46') + 1))
                    continue
                elif ex > 243:
                    egsphant.write(str(Media.index('CORTICAL_BONE_WW86') + 1))
                    continue

            egsphant.write('\n')
        egsphant.write('\n')
    print('...done!')

    print('Writing voxel density assignments...')
    counter = 0
    for z in reversed(range(0, len(ctslice))):
        sys.stdout.write("\r" + str(round(abs(z - m[2]) / m[2] * 100, 2)) + "%")
        sys.stdout.flush()
        for y in range(0, len(ctslice[z].pixel_array)):
            for x in reversed(range(0, len(ctslice[z].pixel_array[y]))):
                ex = ctslice[z].pixel_array[y][x]
                exval = numpy.interp(ex, BigBore_x, BigBore_y)
                if exval < float(0):
                    exval = 0.001225
                if roibladder != -1 and cmap[roibladder][z][x + y * len(xbounds)] == True:
                    exval = 1.00
                if roiurethra != -1 and cmap[roiurethra][z][x + y * len(xbounds)] == True:
                    exval = 1.00
                if counter % 3 == 0:
                    egsphant.write('\n')
                    egsphant.write(' ' * 2)
                egsphant.write("{:.15f}".format(exval))
                if counter % 3 == 0 or counter % 3 == 1:
                    egsphant.write(' ' * 9)
                if counter % 3 == 2:
                    egsphant.write(' ' * 7)
                counter += 1
            egsphant.write('\n')
        egsphant.write('\n')
    print('...done!')

    egsphant.close()

    #############
    # Make EGSINP#
    #############

    egsinp = open(FilePath + 'Patient' + PatientNumber + '_' + labelID + '.egsinp', 'w')
    RAKR = float(rp.SourceSequence[0].ReferenceAirKermaRate)
    AKperHis = float(3.7717 * 10 ** (-14))
    lifetime = float(2056.71)
    DSF = RAKR * lifetime / (AKperHis * 100)

    egsinp.write('Patient' + PatientNumber + ' Prostate egsinp\n\n')
    egsinp.write('#' * 23)
    egsinp.write('\n:Start EGS-mg Geometry:\n\n')
    egsinp.write('1000,1,0\n')
    egsinp.write('6080,' + str(len(seedlocation)) + ',1\n')  # Seed ID should not be hardcoded
    for i in range(0, len(seedlocation)):
        egsinp.write(str(round(
            seedlocation[i].ChannelSequence[0].BrachyControlPointSequence[0].ControlPoint3DPosition[0] * (-1) / 10,
            4)) + ',' + str(
            round(seedlocation[i].ChannelSequence[0].BrachyControlPointSequence[0].ControlPoint3DPosition[1] / 10,
                  4)) + ',' + str(
            round(seedlocation[i].ChannelSequence[0].BrachyControlPointSequence[0].ControlPoint3DPosition[2] / 10,
                  4)) + ',0,0')
        egsinp.write('\n')
    egsinp.write('0,0,0\n\n:Stop EGS-mg Geometry:\n')
    egsinp.write('#' * 23)
    egsinp.write('\n#' + '=' * 79 + '\n')
    egsinp.write('#' * 23)
    egsinp.write(
        '\n:Start Primitive0001:\n\n\nPrimitive Type                  = Rectilinear\nInput X                         '
        '= Individual\nDepth X Boundaries              = -1,1\nInput Y                         = Individual\nDepth Y '
        'Boundaries              = -1,1\nInput Z                         = Individual\nDepth Z Boundaries             '
        ' = -1,1\n\nExtend Outer Dimensions         = On\n\nMedia                           = ' + ','.join(
            Media) + '\n\nSet Medium Everywhere           = 1\n\nDescribe Media By               = Regions\nMedia '
                     'Numbers                   = 1\nFor Media Start Region          = 2\nFor Media Stop Region       '
                     '    = 2\n\nDescribe Important Regions By   = Regions\nFor Important Start Region      = 2\nFor '
                     'Important Stop Region       = 2\n\nIncrease or Decrease Priorities = 1\nFor Priority Start '
                     'Region       = 2\nFor Priority Stop Region        = 2\n\n:Stop Primitive0001:\n')
    egsinp.write('#' * 23)
    egsinp.write('\n#' + '=' * 79 + '\n')
    egsinp.write('#' * 23)
    egsinp.write('\n:Start CT Data Inputs:\n\n')
    egsinp.write(
        'Use CT Data Set to Define Phantom = Yes                                                                      '
        '                             \n')
    egsinp.write(
        'EGSPhant File' + ' ' * 21 + '= ' + '/data/data067/mhaidari/PS/Patient' + PatientNumber + '/Patient' +
        PatientNumber + '_' + labelID + '.egsphant')
    egsinp.write('\nPhantom to be Defined by CT Data  = 1\n\n:Stop CT Data Inputs:\n')
    egsinp.write('#' * 23)
    egsinp.write('\n#' + '=' * 79 + '\n')
    egsinp.write('#' * 23)
    egsinp.write(
        '\n:Start Photon Cross Section Inputs:\nChange Photon Cross Sections           = no\nChange Photon Cross '
        'Sections in Meds   = WATER_0.998\nPercent Change in Photon Cross Section = 0\n\n:Stop Photon Cross Section '
        'Inputs:\n')
    egsinp.write('#' * 23)
    egsinp.write('\n#' + '=' * 79 + '\n')
    egsinp.write('#' * 23)
    egsinp.write(
        '\n:Start I/O Control:\n\nIRestart                     = FIRST\nStore Data Arrays            = Yes\nDose File '
        'Output             = Yes\nOutput Phant File            = NO\nStore Initial Random Numbers = NO\nIWatch       '
        '                = 0\nDose Scaling Factor          = ' + str(
            DSF) + '\n\n:Stop I/O Control:\n')
    egsinp.write('#' * 23)
    egsinp.write('\n#' + '=' * 79 + '\n')
    egsinp.write('#' * 23)
    egsinp.write(
        '\n:Start Monte Carlo Inputs:\n\nInitial Random No. Seeds = 5, 20\nMax CPU Hours Allowed    = 100\nNumber of '
        'Histories	 = 1000000000              \n\n:Stop Monte Carlo Inputs:\n')
    egsinp.write('#' * 23)
    egsinp.write('\n#' + '=' * 79 + '\n')
    egsinp.write('#' * 23)
    egsinp.write(
        '\n:Start MC Transport Parameter:\n\nGlobal ECUT                    = 1.512\nGlobal PCUT                    = '
        '0.001\nGlobal SMAX                    = 1e10\nBound Compton Scattering       = on\nRayleigh Scattering       '
        '     = on\nAtomic Relaxations             = on\nPhotoelectron Angular Sampling = on\nElectron Impact '
        'Ionization     = on   \nBrems Angular Sampling         = KM\nBrems Cross Sections           = nist\nPair '
        'Angular Sampling          = Off   \nESTEPE                         = 0.25\nXIMAX                          = '
        '0.5\nSkin Depth for BCA             = 3.0\nBoundary Crossing Algorithm    = EXACT\nElectron-Step Algorithm   '
        '     = PRESTA-II\nSpin Effects                   = On\nRadiative Compton Corrections  = Off\nPhoton Cross '
        'Sections          = xcom\nFluorescent Photon Cutoff      = 0.001\n\n:Stop MC Transport Parameter:\n')
    egsinp.write('#' * 23)
    egsinp.write('\n#' + '=' * 79 + '\n')
    egsinp.write('#' * 23)
    egsinp.write(
        '\n:Start Muen Data Inputs:\n\nMuendat Input File	 = '
        '/data/data067/mhaidari/egsnrc_MP18/pegs4/data/brachy_xcom_1.5MeV_curr.muendat\n')
    egsinp.write('Require Muendat for Meds = ' + ','.join(Media))
    egsinp.write('\n\n:Stop Muen Data Inputs:\n')
    egsinp.write('#' * 23)
    egsinp.write('\n#' + '=' * 79 + '\n')
    egsinp.write('#' * 23)
    egsinp.write(
        '\n:Start Recycling Inputs:\n\nDo Recycling                      = No\nRotate Recycled Particles         = '
        'No\nTimes to Reuse Recycled Particles = 0\n\n:Stop Recycling Inputs:\n')
    egsinp.write('#' * 23)
    egsinp.write('\n#' + '=' * 79 + '\n')
    egsinp.write('#' * 23)
    egsinp.write(
        '\n:Start Run Mode Inputs:\n\nRun Mode                     = photon source\nPHSP Option                  = NO '
        'PHSP\nPHSP Scoring Object          = 1\nPHSP File Input              =\nScore Dose in Phantoms       = 1\n')
    egsinp.write('Number of non-Source Objects = 1' + '\nObject Containing Sources    = 1\n\n:Stop Run Mode Inputs:\n')
    egsinp.write('#' * 23)
    egsinp.write('\n#' + '=' * 79 + '\n')
    egsinp.write('#' * 23)
    egsinp.write(
        '\n:Start Source Inputs:\n\nIncident Energy              = spectrum\nSpec Filename                = '
        '/home/rthomson/HEN_HOUSE_brachy_distrib/spectra/I125_TG43.spectrum\nSpec IOUTSP                  = '
        'include\nIncident Kinetic Energy(MeV) = 0\nVariable Activity            = Off\nVariable Activity File       '
        '=\n\n:Stop Source Inputs:\n')
    egsinp.write('#' * 23)
    egsinp.write('\n#' + '=' * 79 + '\n')
    egsinp.write('#' * 23)
    egsinp.write(
        '\n:Start Spectrum Scoring Options:\n\nScore Energy Weighted Photon Spectrum Leaving Source = no\n\nScore '
        'Photon Spectrum in a Voxel                     = no\n\nPhantom Containing Spectrum Scoring Voxel            '
        '= 1\nX,Y,Z of Voxel (cm)                                  = 10.0, 0.0, 0.0\nBin Width (MeV)                  '
        '                    = 0.001\n\n:Stop Spectrum Scoring Options:\n')
    egsinp.write('#' * 23)
    egsinp.write('\n#' + '=' * 79 + '\n')
    egsinp.write('#' * 23)
    egsinp.write(
        '\n:Start Volume Correction Inputs:\n\nVolume Correction                                         = correct '
        'volume\nCorrect Voxel Volumes in Phantom                          = 1\nDensity of Random Points for Vol Cor '
        '(cm^-3)              = 1E7\nDensity of Random Points for Finding Seed Regions (cm^-3) = 1.E6\nAdditional '
        'Volume Correction Needed                       = no\nRegion (XMin,XMax,YMin,YMax,ZMin,ZMax)                  '
        '  = -1.0,1.0,-1.0,1.0,-0.3,0.5\n\n:Stop Volume Correction Inputs:\n')
    egsinp.write('#' * 23)

    egsinp.close()
