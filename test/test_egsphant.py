#!/usr/bin/env python3
import unittest
from egsphant.egsphant import EGSphant


class TestEGSphant(unittest.TestCase):

    def setUp(self):
        self.egsphant = EGSphant('test_destination.egsphant')

    def test_number_of_media(self):
        self.assertEqual(self.egsphant.number_of_media, 6, 'incorrect number of media')

    def test_media(self):
        media = ['Medium_1', 'Medium_2', 'Medium_3', 'Medium_4', 'Medium_5', 'Medium_6']
        self.assertListEqual(self.egsphant.media, media, 'incorrect list of media')

    def test_medium_keys(self):
        medium_keys = {'1': 'Medium_1', '2': 'Medium_2', '3': 'Medium_3', '4': 'Medium_4', '5': 'Medium_5',
                       '6': 'Medium_6'}
        self.assertDictEqual(self.egsphant.medium_keys, medium_keys, 'incorrect dictionary of medium keys')

    def test_dimensions(self):
        self.assertListEqual(self.egsphant.dimensions, [3, 4, 5], 'incorrect dimensions')

    def test_bounds(self):
        bounds = [[-3, -1, 1, 3],
                  [-2, -1, 0, 1, 2],
                  [-1.5, -1, -0.5, 0.5, 1, 1.5]]
        self.assertListEqual(self.egsphant.bounds, bounds, 'incorrect bounds')

    def test_phantom_array(self):
        phantom_flat = [1, 1, 1, 1, 2, 1, 1, 2, 1, 1, 1, 1,
                        3, 3, 3, 3, 2, 3, 3, 2, 3, 3, 3, 3,
                        4, 4, 4, 4, 5, 4, 4, 5, 4, 4, 4, 4,
                        6, 6, 6, 6, 5, 6, 6, 5, 6, 6, 6, 6,
                        1, 2, 3, 4, 5, 6, 3, 2, 1, 6, 5, 4]
        phantom = list(map(str, phantom_flat))
        self.assertListEqual(self.egsphant.phantom.flatten(order='F').tolist(), phantom, 'incorrect phantom array')

    def test_density_array(self):
        density = [float(i) for i in range(60)]
        self.assertListEqual(self.egsphant.density.flatten(order='F').tolist(), density, 'incorrect density array')

    def test_nvox(self):
        self.assertEqual(self.egsphant.nvox(), 60, 'incorrect number of voxels')

    def test_get_medium_key(self):
        self.assertEqual(self.egsphant.get_medium_key('Medium_2'), '2', 'incorrect medium key')

    def test_get_index_voxels_with_medium(self):
        self.assertListEqual(self.egsphant.get_index_of_voxels_with_medium('Medium_1'),
                             [0, 1, 2, 3, 5, 6, 8, 9, 10, 11, 48, 56],
                             'incorrect index of voxels for Medium_1')

    # def test_write_egsphant(self):
    #     self.egsphant.write_egsphant('write.egsphant')
    #     with open('test_destination.egsphant', 'r') as file:
    #         expected_lines = [line.split() for line in file]
    #     with open('write.egsphant', 'r') as file:
    #         written_lines = [line.split() for line in file]
    #     self.assertListEqual(expected_lines, written_lines, 'did not write to file correctly')

    def test_create_water_egsphant(self):
        bounds = [[-3, -2, -1, 1, 2, 3],
                  [-1.5, -1, -0.5, 0, 0.5, 1, 1.5],
                  [-2, 1, 0, 1, 2]]
        egsphant = EGSphant.create_water_egsphant(bounds)
        self.assertEqual(egsphant.number_of_media, 1, 'incorrect number of media for water phantom')
        self.assertListEqual(egsphant.media, ['WATER_0.998'], 'incorrect water medium')
        self.assertDictEqual(egsphant.medium_keys, {'1': 'WATER_0.998'}, 'incorrect medium key for water')
        self.assertListEqual(egsphant.dimensions, [5, 6, 4], 'incorrect water phantom dimensions')
        self.assertListEqual(egsphant.bounds, bounds, 'incorrect water phantom bounds')

if __name__ == '__main__':
    unittest.main()
