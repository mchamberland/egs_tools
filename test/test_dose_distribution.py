#!/usr/bin/env python3
import unittest
from dose_distribution.dose3d import DoseDistribution


class TestDoseDistribution(unittest.TestCase):

    def setUp(self):
        self.dose = DoseDistribution('read.3ddose')

        with open('read.3ddose', 'r') as file:
            lines = [line.split() for line in file]
            self.expected_dose = list(map(float, lines[4]))
            self.expected_fract_unc = list(map(float, lines[5]))

    def test_read_dimensions(self):
        rect_dimensions = [3, 2, 5]
        self.assertEqual(self.dose.dimensions, rect_dimensions,
                         'incorrect dimensions list')

    def test_read_x_bounds(self):
        x_bounds = [-3.0, -1.0, 1.0, 3.0]
        self.assertEqual(self.dose.bounds[0], x_bounds,
                         'incorrect x_bounds')

    def test_read_y_bounds(self):
        y_bounds = [-1.0, 0.0, 1.0]
        self.assertEqual(self.dose.bounds[1], y_bounds,
                         'incorrect y_bounds')

    def test_read_z_bounds(self):
        z_bounds = [-10.0, -6.0, -2.0, 2.0, 6.0, 10.0]
        self.assertEqual(self.dose.bounds[2], z_bounds,
                         'incorrect z_bounds')

    def test_read_dose_array(self):
        dose = self.dose.dose.flatten('F').tolist()
        self.assertListEqual(dose, self.expected_dose,
                             'incorrect nominal dose values')

    def test_read_unc_array(self):
        fract_unc = self.dose.fract_unc.flatten('F').tolist()
        for (u, exp) in zip(fract_unc, self.expected_fract_unc):
            self.assertAlmostEqual(u, exp, 7,
                                   'incorrect uncertainty values')

    def test_write_dose(self):
        self.dose.write_dose('write.3ddose')
        with open('read.3ddose', 'r') as file:
            expected_lines = [line.split() for line in file]
        with open('write.3ddose', 'r') as file:
            written_lines = [line.split() for line in file]
        self.assertListEqual(expected_lines, written_lines, 'did not write to file correctly')

    def test_nvox(self):
        self.assertEqual(self.dose.nvox(), 30,
                         'incorrect number of voxels for distribution')

    def test_get_max_dose(self):
        max_dose, unc = 29.66447663, 0.00156969
        self.assertAlmostEqual(max_dose, self.dose.get_max_dose()[0], 7, 'incorrect max dose')
        self.assertAlmostEqual(unc, self.dose.get_max_dose()[1], 7, 'incorrect unc on max dose')

    def test_get_min_dose(self):
        min_dose, unc = 1.19834069, 0.00043182
        self.assertAlmostEqual(min_dose, self.dose.get_min_dose()[0], 7, 'incorrect min dose')
        self.assertAlmostEqual(unc, self.dose.get_min_dose()[1], 7, 'incorrect unc on min dose')


if __name__ == '__main__':
    unittest.main()
