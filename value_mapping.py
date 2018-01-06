# This algorithm taken from Motion Imagery Standards Board, MISB ST 1201.1
# Standard: Floating Point to Integer Mapping
# http://www.gwg.nga.mil/misb/docs/standards/ST1201.1.pdf
import math


class FloatingPointToIntegerMapping:
    def __init__(self, min_float, max_float, precision_or_bytes, mode='bytes'):
        self.forward_mapping_factor = None
        self.reverse_mapping_factor = None
        self.min = min_float
        self.max = max_float

        if mode == 'precision':
            bits = math.ceil(math.log((max_float - min_float) / precision_or_bytes, 2)) + 1
            num_bytes = math.ceil(bits / 8)
        elif mode == 'bytes':
            num_bytes = precision_or_bytes
        else:
            raise Exception('Unknown mode. Use ''bytes'' or ''precision''.')

        bpow = math.ceil(math.log(max_float - min_float, 2))
        dpow = 8 * num_bytes - 1
        self.forward_mapping_factor = math.pow(2, dpow - bpow)
        self.reverse_mapping_factor = math.pow(2, bpow - dpow)
        self.offset = 0

        if min_float < 0 < max_float:
            self.offset = self.forward_mapping_factor * min_float - math.floor(self.forward_mapping_factor * min_float)

    def float_to_integer(self, x):
        return math.floor(self.forward_mapping_factor * (x - self.min) + self.offset)

    def integer_to_float(self, x):
        return self.reverse_mapping_factor * (x - self.offset) + self.min
