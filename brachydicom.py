from os import listdir
from os.path import join, isfile, splitext
import pydicom


def get_all_files_in_directory(directory='.'):
    return [f for f in listdir(directory) if isfile(join(directory, f))]


def get_dicom_files_in_directory(directory='.'):
    all_files = get_all_files_in_directory(directory)
    return [f for f in all_files if f.endswith('.dcm')]


def get_plan_files_in_directory(directory='.'):
    dicom_files = get_dicom_files_in_directory(directory)
    return [f for f in dicom_files if f.startswith('RP')]


def get_structure_files_in_directory(directory='.'):
    dicom_files = get_dicom_files_in_directory(directory)
    return [f for f in dicom_files if f.startswith('RS')]


def get_dose_files_in_directory(directory='.'):
    dicom_files = get_dicom_files_in_directory(directory)
    return [f for f in dicom_files if f.startswith('RD')]


def open_dicom_files(dicom_files):
    filenames = [splitext(file)[0] for file in dicom_files]
    dicom_files = [pydicom.read_file(file) for file in dicom_files]
    return filenames, dicom_files


def get_dwells(dicom_file):
    root = dicom_file.ApplicationSetupSequence[0].ChannelSequence
    number_of_channels = len(root)
    dwell_positions = []
    dwell_times = []

    for i in range(number_of_channels):
        number_of_control_points = int(root[i].NumberOfControlPoints)
        channel_total_time = float(root[i].ChannelTotalTime)
        final_cumulative_time_weight = float(root[i].FinalCumulativeTimeWeight)
        ctrl_pts = root[i].BrachyControlPointSequence

        current_time = 0
        current_position = [-9999, -9999, -9999]

        for j in range(number_of_control_points):
            previous_position = current_position
            current_position = list(map(float, ctrl_pts[j].ControlPoint3DPosition))

            previous_time = current_time
            cumulative_time_weight = float(ctrl_pts[j].CumulativeTimeWeight)
            current_time = channel_total_time * cumulative_time_weight / final_cumulative_time_weight

            if current_position == previous_position:
                dwell_positions.append(current_position)
                dwell_times.append(current_time - previous_time)

    return dwell_positions, dwell_times, number_of_channels
