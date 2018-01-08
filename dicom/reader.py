from os import listdir
from os.path import join, isfile, splitext
import pydicom
import numpy as np
from typing import Tuple
import datetime

SOURCE_MODELS = {'I-125 (6711) [TG43 updated]': 'OncoSeed_6711',
                 'I-125(AgX100)': 'AgX100'}


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


def read_ct_files_in_directory(directory='.'):
    files = list_ct_files_in_directory(directory)
    filenames = [splitext(file)[0] for file in files]
    ct_files = [pydicom.read_file(join(directory, file)) for file in files]
    return filenames, ct_files


def get_sorted_ct_files_in_directory(directory='.'):
    # sorted by z-slice
    filenames, ct_files = read_ct_files_in_directory(directory)
    zslices = [float(ct.ImagePositionPatient[2]) for ct in ct_files]
    temp = list(zip(zslices, filenames, ct_files))
    temp.sort()
    zslices, filenames, ct_files = zip(*temp)
    return filenames, ct_files


def read_plan_files_in_directory(directory='.'):
    files = list_plan_files_in_directory(directory)
    filenames = [splitext(file)[0] for file in files]
    dicom_plans = [pydicom.read_file(join(directory, file)) for file in files]
    return filenames, dicom_plans


def read_dose_files_in_directory(directory='.'):
    files = list_dose_files_in_directory(directory)
    filenames = [splitext(file)[0] for file in files]
    dicom_doses = [pydicom.read_file(join(directory, file)) for file in files]
    return filenames, dicom_doses


def read_structure_files_in_directory(directory='.'):
    files = list_structure_files_in_directory(directory)
    filenames = [splitext(file)[0] for file in files]
    dicom_structs = [pydicom.read_file(join(directory, file)) for file in files]
    return filenames, dicom_structs


def get_all_dwells_in_mm_and_s(dicom_file) -> Tuple[np.ndarray, np.ndarray, int]:
    # these include dwells that have no source, to use possibly to determine source orientation in 3D
    # positions in mm and times in s
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

    return np.array(dwell_positions), np.array(dwell_times), number_of_channels


def get_dwells_in_mm_and_s(dicom_file) -> Tuple[np.ndarray, np.ndarray, int]:
    # positions in mm and times in s
    dwells, times, number_of_channels = get_all_dwells_in_mm_and_s(dicom_file)
    dwell_times = times[times > 0]
    dwell_positions = dwells[times > 0]
    return dwell_positions, dwell_times, number_of_channels


def get_modality(dicom_file):
    return dicom_file.Modality


def get_treatment_type(dicom_file):
    return dicom_file.BrachyTreatmentType


def get_treatment_technique(dicom_file):
    return dicom_file.BrachyTreatmentTechnique


def get_rakr_in_ugy_per_h_at_1m(dicom_file):
    # in microGy per hour at 1 m
    return float(dicom_file.SourceSequence[0].ReferenceAirKermaRate)


def get_half_life_in_days(dicom_file):
    # in days
    return float(dicom_file.SourceSequence[0].SourceIsotopeHalfLife)


def get_source_strength_reference_date_and_time(dicom_file):
    date = dicom_file.SourceSequence[0].SourceStrengthReferenceDate
    time = dicom_file.SourceSequence[0].SourceStrengthReferenceTime
    date_time = date + time
    return datetime.datetime.strptime(date_time, "%Y%m%d%H%M%S")


def get_study_date_and_time(dicom_file):
    date = dicom_file.StudyDate
    time = dicom_file.StudyTime
    date_time = date + time
    if not time:
        return datetime.datetime.strptime(date_time, "%Y%m%d")
    else:
        return datetime.datetime.strptime(date_time, "%Y%m%d%H%M%S")


def get_total_treatment_time_in_s(dicom_file):
    # in seconds
    rakr = get_rakr_in_ugy_per_h_at_1m(dicom_file)
    total_reference_air_kerma = float(dicom_file.SourceSequence[0].ApplicationSetupSequence[0].TotalReferenceAirKerma)
    return total_reference_air_kerma / rakr * 3600


def get_source_model(dicom_file):
    return SOURCE_MODELS[dicom_file.SourceSequence[0].SourceIsotopeName]


def get_number_of_seeds(dicom_file):
    return len(dicom_file.ApplicationSetupSequence)


def get_seed_locations_in_mm(dicom_file) -> np.ndarray:
    # positions in mm
    seed_positions = []

    for seed in dicom_file.ApplicationSetupSequence:
        pos = [float(x) for x in seed.ChannelSequence[0].BrachyControlPointSequence[0].ControlPoint3DPosition]
        seed_positions.append(pos)

    return np.array(seed_positions)
