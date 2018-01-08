import sys


class RTDoseInfo:
    def __init__(self, dicom):
        if dicom.Modality != 'RTDOSE':
            raise Exception('This class only accepts DICOM RT Dose objects.')
        self.dicom = dicom
        self.dx, self.dy = self.get_pixel_size()
        self.dz = float(dicom.GridFrameOffsetVector[1])
        self.z_direction = 1 if self.dz > 0 else -1
        self.nx, self.ny, self.nz = self.get_dimensions()
        self.centre_of_first_pixel = self.get_centre_of_first_pixel()
        self.first_slice_position = float(dicom.ImagePositionPatient[2])
        self.grid_extents = self.get_grid_extents()
        self.grid_size = self.get_grid_size()
        self.centre = self.get_centre()
        self.dose_max = self.get_max_dose()
        self.dose_units = self.get_dose_units()

    def get_pixel_size(self):
        return float(self.dicom.PixelSpacing[0]), float(self.dicom.PixelSpacing[1])

    def get_dimensions(self):
        return self.dicom.Columns, self.dicom.Rows, int(self.dicom.NumberOfFrames)

    def get_centre_of_first_pixel(self):
        return float(self.dicom.ImagePositionPatient[0]), float(self.dicom.ImagePositionPatient[1])

    def get_grid_extents(self):
        x_initial = self.centre_of_first_pixel[0] - self.dx / 2
        y_initial = self.centre_of_first_pixel[1] - self.dy / 2
        z_initial = self.first_slice_position - self.dz / 2

        x_final = x_initial + self.nx * self.dx
        y_final = y_initial + self.ny * self.dy
        z_final = z_initial + self.nz * self.dz

        return (x_initial, x_final), (y_initial, y_final), (z_initial, z_final)

    def get_grid_size(self):
        (x_initial, x_final), (y_initial, y_final), (z_initial, z_final) = self.grid_extents

        return x_final - x_initial, y_final - y_initial, abs(z_final - z_initial)

    def get_centre(self):
        return (self.grid_extents[0][0] + self.grid_size[0] / 2,
                self.grid_extents[1][0] + self.grid_size[1] / 2,
                self.grid_extents[2][0] + self.z_direction * self.grid_size[2] / 2)

    def get_max_dose(self):
        return (self.dicom.pixel_array * self.dicom.DoseGridScaling).max()

    def get_dose_units(self):
        return self.dicom.DoseUnits.lower().capitalize()

    def print_info(self, save_to_file=None):
        if save_to_file:
            file = open(save_to_file + '.rtdose_info', 'w')
        else:
            file = sys.stdout
        print("Dose grid dimensions (voxels):\n{} x {} x {}\n".format(self.nx, self.ny, self.nz), file=file)
        print("Voxel size:\n{} mm x {} mm x {} mm\n".format(self.dx, self.dy, abs(self.dz)), file=file)
        print("Extents of dose grid:", file=file)
        print("{:.3f} cm to {:.3f} cm along x.".format(self.grid_extents[0][0] / 10, self.grid_extents[0][1] / 10),
              file=file)
        print("{:.3f} cm to {:.3f} cm along y.".format(self.grid_extents[1][0] / 10, self.grid_extents[1][1] / 10),
              file=file)
        print("{:.3f} cm to {:.3f} cm along z.\n".format(self.grid_extents[2][0] / 10, self.grid_extents[2][1] / 10),
              file=file)
        print("Centre of dose grid:\n({:.2f} cm, {:.2f} cm, {:.2f} cm)\n".format(self.centre[0] / 10,
                                                                                 self.centre[1] / 10,
                                                                                 self.centre[2] / 10), file=file)
        print("Total size of dose grid:\n{:.2f} cm x {:.2f} cmx {:.2f} cm\n".format(self.grid_size[0] / 10,
                                                                                    self.grid_size[1] / 10,
                                                                                    self.grid_size[2] / 10), file=file)

        print("Max dose:\n{:.2f} {}".format(self.dose_max, self.dose_units), file=file)
