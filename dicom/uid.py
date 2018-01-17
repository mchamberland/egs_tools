import os
import datetime as dt

DICOM_SOP_UID = {'CT': '1.2.840.10008.5.1.4.1.1.2',
                 'RTDOSE': '1.2.840.10008.5.1.4.1.1.481.2',
                 'RTSTRUCT': '1.2.840.10008.5.1.4.1.1.481.3',
                 'RTPLAN': '1.2.840.10008.5.1.4.1.1.481.5'}

GENERIC_ROOT_UID = "1.2.124"


class DicomUIDGenerator:
    def __init__(self, dicom_type=None):
        if dicom_type:
            dicom_type = dicom_type.upper()
        if dicom_type in DICOM_SOP_UID:
            self.sop = DICOM_SOP_UID[dicom_type]
        else:
            self.sop = None
        self.root = GENERIC_ROOT_UID
        self.uid_counter = 0

    def create_series_instance_uid(self):
        return self.create_uid(self.root + '.1.3')

    def create_sop_instance_uid(self):
        if self.sop:
            return self.sop
        else:
            return self.create_uid(self.root + '.1.4')

    def create_study_instance_uid(self):
        return self.create_uid(self.root + '.1.2')

    def create_frame_of_reference_uid(self):
        return self.create_uid(self.root + '.1.5')

    def create_uid(self, uid):
        uid = uid + '.' + str(os.getpid())
        now = dt.datetime.strftime(dt.datetime.today(), "%Y%m%d%H%M%f")
        uid = uid + '.' + now
        self.uid_counter += 1
        return uid + '.' + str(self.uid_counter)


def generate_new_uids(dicom_type=None):
    generator = DicomUIDGenerator(dicom_type)
    sop = generator.create_sop_instance_uid()
    study = generator.create_study_instance_uid()
    series = generator.create_series_instance_uid()
    frame = generator.create_frame_of_reference_uid()

    return sop, study, series, frame
