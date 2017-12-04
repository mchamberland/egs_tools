import os


class CTConversionToTissue:
    def __init__(self, directory="."):
        self.directory = directory
        self.media_list = []

    def read_ctconv_file(self, filename='default.ctconv'):
        if not os.path.exists(filename):
            return -1

        with open(filename, 'r') as file:
            lines = file.readlines()
            lines = [line.strip() for line in lines if line.strip()]

        for line in lines:
            name, density, min_ctnum, max_ctnum = [x for x in line]
            medium = MediumInfo(name, float(density), int(min_ctnum), int(max_ctnum))
            self.media_list.append(medium)


class MediumInfo:
    def __init__(self, name=None, density=None, min_ctnum=-1000, max_ctnum=3000):
        self.name = name
        self.density = density
        self.min_ctnum = min_ctnum
        self.max_ctnum = max_ctnum
