# Testing code

import numpy as np
import unittest
from .. import insar_vector_functions


class Tests(unittest.TestCase):

    def test_lkv_conversion(self):
        """ Test that look vectors are converted to azimuth/incidence, for scalar and numpy array operation"""
        [lkv_e, lkv_n, lkv_u] = insar_vector_functions.flight_incidence_angles2look_vector(190, 30)
        [out_flight, out_inc] = insar_vector_functions.look_vector2flight_incidence_angles(lkv_e, lkv_n, lkv_u)
        self.assertAlmostEqual(190, out_flight)
        self.assertAlmostEqual(30, out_inc)

        incidence_grd = np.multiply(37, np.ones((10, 10)))  # in degrees
        azimuth_grd = np.multiply(190, np.ones((10, 10)))  # in degrees
        e_grd, n_grd, u_grd = insar_vector_functions.flight_incidence_angles2look_vector(azimuth_grd, incidence_grd)
        azimuth_ret, incidence_ret = insar_vector_functions.look_vector2flight_incidence_angles(e_grd, n_grd, u_grd)
        np.testing.assert_allclose(azimuth_grd, azimuth_ret)
        np.testing.assert_allclose(incidence_grd, incidence_ret)
        return

    def test_proj_into_vertical_no_horizontal(self):
        # LKV 45 degree test case:
        semivert = insar_vector_functions.proj_los_into_vertical_no_horiz(1, [0.7, 0, 0.7])
        self.assertAlmostEqual(semivert, 1.41421356)
        # LKV vertical test case:
        semivert = insar_vector_functions.proj_los_into_vertical_no_horiz(1, [0, 0, 1])
        self.assertAlmostEqual(semivert, 1.0)
        return


if __name__ == "__main__":
    unittest.main()
