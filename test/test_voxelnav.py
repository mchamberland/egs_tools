#!/usr/bin/env python3
import unittest
import numpy as np
import voxelnav


class TestVoxelnav(unittest.TestCase):

    def setUp(self):
        self.bounds = [[-3, -1, 1, 3],
                       [-2, -1, 0, 1, 2],
                       [-1.5, -1, -0.5, 0.5, 1, 1.5]]

        self.voxels_flat = [1, 1, 1, 1, 2, 1, 1, 2, 1, 1, 1, 1,
                            3, 3, 3, 3, 2, 3, 3, 2, 3, 3, 3, 3,
                            4, 4, 4, 4, 5, 4, 4, 5, 4, 4, 4, 4,
                            6, 6, 6, 6, 5, 6, 6, 5, 6, 6, 6, 6,
                            1, 2, 3, 4, 5, 6, 3, 2, 1, 6, 5, 4]

        self.dimensions = [3, 4, 5]

        voxels = np.array(self.voxels_flat)
        self.voxels = voxels.reshape(self.dimensions, order='F')

    def test_get_voxel_center_from_index(self):
        self.assertTupleEqual(voxelnav.get_voxel_center_from_index(13, self.bounds),
                              (0, -1.5, -0.75), 'incorrect voxel centre from index')

    def test_get_voxel_center_from_ijk(self):
        self.assertTupleEqual(voxelnav.get_voxel_center_from_ijk((2, 2, 2), self.bounds),
                              (2, 0.5, 0), 'incorrect voxel centre from multi-index')

    def test_get_voxel_size_from_ijk(self):
        self.assertEqual(voxelnav.get_voxel_size_from_ijk((2, 2, 2), self.bounds),
                         (2, 1, 1), 'incorrect size from ijk')

    def test_get_index_from_ijk(self):
        self.assertEqual(voxelnav.get_index_from_ijk((2, 2, 2), self.dimensions),
                         32, 'incorrect index from ijk')

    def test_get_ijk_from_index(self):
        self.assertTupleEqual(voxelnav.get_ijk_from_index(32, self.dimensions),
                              (2, 2, 2), 'incorrect ijk from index')

    def test_get_index_from_position(self):
        self.assertEqual(voxelnav.get_index_from_position(0.2, self.bounds[0]),
                         1, 'incorrect index from 1-D position')
        self.assertEqual(voxelnav.get_index_from_position(-3.5, self.bounds[0], is_strictly_inside=False),
                         0, 'incorrect index from 1-D position')
        self.assertEqual(voxelnav.get_index_from_position(3.5, self.bounds[0], is_strictly_inside=False),
                         2, 'incorrect index from 1-D position')

    def test_region_indices_from_low_and_hi_pos(self):
        self.assertTupleEqual(voxelnav.get_region_indices_from_low_and_hi_pos(-4, -3.5, self.bounds[0]),
                              (-1, -1), 'incorrect region indices from outside bounds')
        self.assertTupleEqual(voxelnav.get_region_indices_from_low_and_hi_pos(-2, 2, self.bounds[0]),
                              (0, 2), 'incorrect region indices from positions')
        self.assertTupleEqual(voxelnav.get_region_indices_from_low_and_hi_pos(3, 3.5, self.bounds[0]),
                              (2, 2), 'incorrect region indices from edge case')
        self.assertTupleEqual(voxelnav.get_region_indices_from_low_and_hi_pos(-3, 3.5, self.bounds[0]),
                              (0, 2), 'incorrect region indices from edge case')
        self.assertTupleEqual(voxelnav.get_region_indices_from_low_and_hi_pos(-3, 3, self.bounds[0]),
                              (0, 2), 'incorrect region indices from edge case')
        self.assertTupleEqual(voxelnav.get_region_indices_from_low_and_hi_pos(0, 3, self.bounds[0]),
                              (1, 2), 'incorrect region indices from edge case')
        self.assertTupleEqual(voxelnav.get_region_indices_from_low_and_hi_pos(0, 3, self.bounds[0],
                                                                              include_partial_voxel=False),
                              (2, 2), 'incorrect region indices from edge case')
        self.assertTupleEqual(voxelnav.get_region_indices_from_low_and_hi_pos(-2, 2, self.bounds[0],
                                                                              include_partial_voxel=False),
                              (1, 1), 'incorrect region indices from edge case')

    def test_get_next_x_neighbour(self):
        self.assertEqual(voxelnav.get_next_x_neighbour((0, 1, 2), self.bounds[0]),
                         1, 'incorrect next x neighbour')
        self.assertEqual(voxelnav.get_next_x_neighbour((2, 1, 2), self.bounds[0]),
                         2, 'incorrect next x neighbour')
        self.assertEqual(voxelnav.get_next_x_neighbour((3, 1, 2), self.bounds[0]),
                         -1, 'incorrect next x neighbour')

    def test_get_prev_x_neighbour(self):
        self.assertEqual(voxelnav.get_prev_x_neighbour((0, 1, 2)),
                         0, 'incorrect prev x neighbour')
        self.assertEqual(voxelnav.get_prev_x_neighbour((2, 1, 2)),
                         1, 'incorrect prev x neighbour')

    def test_get_next_y_neighbour(self):
        self.assertEqual(voxelnav.get_next_y_neighbour((0, 1, 2), self.bounds[1]),
                         2, 'incorrect next y neighbour')
        self.assertEqual(voxelnav.get_next_y_neighbour((2, 3, 2), self.bounds[1]),
                         3, 'incorrect next y neighbour')
        self.assertEqual(voxelnav.get_next_y_neighbour((3, 4, 2), self.bounds[1]),
                         -1, 'incorrect next y neighbour')

    def test_get_prev_y_neighbour(self):
        self.assertEqual(voxelnav.get_prev_y_neighbour((0, 0, 2)),
                         0, 'incorrect prev y neighbour')
        self.assertEqual(voxelnav.get_prev_y_neighbour((2, 1, 2)),
                         0, 'incorrect prev y neighbour')

    def test_get_next_z_neighbour(self):
        self.assertEqual(voxelnav.get_next_z_neighbour((0, 1, 2), self.bounds[2]),
                         3, 'incorrect next z neighbour')
        self.assertEqual(voxelnav.get_next_z_neighbour((2, 3, 4), self.bounds[2]),
                         4, 'incorrect next z neighbour')
        self.assertEqual(voxelnav.get_next_z_neighbour((3, 4, 5), self.bounds[2]),
                         -1, 'incorrect next z neighbour')

    def test_get_prev_z_neighbour(self):
        self.assertEqual(voxelnav.get_prev_z_neighbour((0, 0, 0)),
                         0, 'incorrect prev z neighbour')
        self.assertEqual(voxelnav.get_prev_z_neighbour((2, 1, 2)),
                         1, 'incorrect prev z neighbour')

if __name__ == '__main__':
    unittest.main()
