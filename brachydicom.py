from os import listdir
from os.path import join, isfile, splitext
import pydicom


def get_dwell_positions_and_times():

    all_files = [f for f in listdir('.') if isfile(join('.', f))]
    dicom_files = [f for f in all_files if f.endswith('.dcm')]
    plan_files = [f for f in dicom_files if f.startswith('RP')]

    for file in plan_files:
        dc = pydicom.read_file(file)
        root = dc.ApplicationSetupSequence[0].ChannelSequence
        length = len(root)
        dwell_filename = splitext(file)[0] + '.dwell'

        with open(dwell_filename, 'wt') as dwell_file:
            dwell_pos_and_times = ''

            # TODO should take all dwell positions to determine direction cosines and, hence, rotation of source!
            for i in range(length):
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

                    # if the positions are different or the times are identical, then the source is in transit
                    if current_position != previous_position or current_time == previous_time:
                        pass
                    # this is a dwell position
                    else:
                        time = current_time - previous_time
                        x, y, z = current_position
                        dwell_pos_and_times +=\
                            '{0:.2f}\t{1:.2f}\t{2:.2f}\t{3:.2f}\n'.format(x / 10, y / 10, z / 10, time)

            dwell_file.write(dwell_pos_and_times)
