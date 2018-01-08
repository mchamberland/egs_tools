import os
import datetime as dt

ROOT_UID = "1.2.124"


class DicomUIDGenerator:
    def __init__(self, root=None):
        if root:
            self.root = root
        else:
            self.root = ROOT_UID
        self.uid_counter = 0

    def create_series_instance_uid(self):
        return self.create_uid(self.root + '.1.3')

    def create_sop_instance_uid(self):
        return self.create_uid(self.root + '1.4')

    def create_study_instance_uid(self):
        return self.create_uid(self.root + '1.2')

    def create_frame_of_reference_uid(self):
        return self.create_uid(self.root + '1.5')

    def create_uid(self, uid):
        uid = uid + '.' + str(os.getpid())
        now = dt.datetime.strftime(dt.datetime.today(), "%Y%m%d%H%M%f")
        uid = uid + '.' + now
        self.uid_counter += 1
        return uid + '.' + str(self.uid_counter)
