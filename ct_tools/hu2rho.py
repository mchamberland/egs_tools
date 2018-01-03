import os
from os.path import join
from scipy.interpolate import interp1d

PATH_TO_CT_CALIBRATION = os.path.expandvars("$CT_CALIBRATION_DIR")


class HU2rho:
    def __init__(self, filename, directory=PATH_TO_CT_CALIBRATION):
        self.ct_numbers = []
        self.densities = []
        path = join(directory, filename)
        if os.path.exists(path):
            self.load_hu_to_density_calibration(path)
        else:
            raise FileNotFoundError

    def load_hu_to_density_calibration(self, filename):
        with open(filename, 'r') as file:
            lines = file.readlines()
            ctnum_rho_pairs = [line.split() for line in lines]
            for ctnum, rho in ctnum_rho_pairs:
                self.ct_numbers.append(float(ctnum))
                self.densities.append(float(rho))

    def get_density_from_hu(self, ctnum, extrapolate=False):
        f = interp1d(self.ct_numbers, self.densities, fill_value='extrapolate')
        if extrapolate:
            density = float(f(ctnum))
        else:
            if ctnum < self.ct_numbers[0]:
                ctnum = self.ct_numbers[0]
            elif ctnum > self.ct_numbers[-1]:
                ctnum = self.ct_numbers[-1]
            density = float(f(ctnum))
        if density > 0:
            return density
        else:
            return 0
