#!/usr/bin/env python3
import unittest
from egsphant.egsphant import EGSphant
import egsphant.manip


class TestEgsphantManip(unittest.TestCase):

    def setUp(self):
        self.egsphant = EGSphant('test_destination.egsphant')

    def test_trim_phantom(self):
        egsphant.manip.trim_phantom(self.egsphant, (-3.5, 3.5, -0.5, 0.5, -1.5, -0.5))
        self.egsphant.write_egsphant('trimmed.egsphant')

        phantom_flat = ['1', '2', '1', '1', '2', '1', '3', '2', '3', '3', '2', '3']
        density_flat = [3, 4, 5, 6, 7, 8, 15, 16, 17, 18, 19, 20]
        self.assertListEqual(self.egsphant.dimensions, [3, 2, 2], 'incorrect trimmed dimensions')
        self.assertListEqual(self.egsphant.bounds[1], [-1, 0, 1], 'incorrect y-bounds')
        self.assertListEqual(self.egsphant.bounds[2], [-1.5, -1, -0.5], 'incorrect z-bounds')
        self.assertListEqual(self.egsphant.phantom.flatten(order='F').tolist(), phantom_flat,
                             'incorrect trimmed phantom')
        self.assertListEqual(self.egsphant.density.flatten(order='F').tolist(), density_flat,
                             'incorrect trimmed density')

    def test_copy_medium_from_source_to_destination(self):
        source = EGSphant('test_source.egsphant')
        egsphant.manip.copy_medium_from_source_to_destination('Medium2_5', source, self.egsphant)
        self.egsphant.write_egsphant('copied_Medium2_5.egsphant')
        # TODO test egsphants of different dimensions/bounds

    def test_replace_original_medium_with_new_medium(self):
        egsphant.manip.replace_original_medium_with_new_medium(self.egsphant, 'Medium_1', 'New_medium', 0.99)
        self.egsphant.write_egsphant('replace_Medium_1_by_New_medium.egsphant')

    def test_add_medium(self):
        pass
        # TODO test_add_medium

if __name__ == '__main__':
    unittest.main()
