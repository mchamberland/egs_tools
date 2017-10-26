from os import listdir
from os.path import join, isfile, splitext
import pydicom


def list_all_files_in_directory(directory='.'):
    return [f for f in listdir(directory) if isfile(join(directory, f))]


def list_dicom_files_in_directory(directory='.'):
    all_files = list_all_files_in_directory(directory)
    return [f for f in all_files if f.endswith('.dcm')]


def list_plan_files_in_directory(directory='.'):
    dicom_files = list_dicom_files_in_directory(directory)
    return [f for f in dicom_files if pydicom.read_file(join(directory, f)).Modality == 'RTPLAN']


def list_structure_files_in_directory(directory='.'):
    dicom_files = list_dicom_files_in_directory(directory)
    return [f for f in dicom_files if pydicom.read_file(join(directory, f)).Modality == 'RTSTRUCT']


def list_dose_files_in_directory(directory='.'):
    dicom_files = list_dicom_files_in_directory(directory)
    return [f for f in dicom_files if pydicom.read_file(join(directory, f)).Modality == 'RTDOSE']


def list_ct_files_in_directory(directory='.'):
    dicom_files = list_dicom_files_in_directory(directory)
    return [f for f in dicom_files if pydicom.read_file(join(directory, f)).Modality == 'CT']


def read_plan_files_in_directory():
    files = list_plan_files_in_directory()
    filenames = [splitext(file)[0] for file in files]
    dicom_plans = [pydicom.read_file(file) for file in files]
    return filenames, dicom_plans


def read_dose_files_in_directory():
    files = list_dose_files_in_directory()
    filenames = [splitext(file)[0] for file in files]
    dicom_doses = [pydicom.read_file(file) for file in files]
    return filenames, dicom_doses


def read_structure_files_in_directory():
    files = list_structure_files_in_directory()
    filenames = [splitext(file)[0] for file in files]
    dicom_structs = [pydicom.read_file(file) for file in files]
    return filenames, dicom_structs


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


def read_modality(dicom_file):
    modality = dicom_file.Modality

    return modality
