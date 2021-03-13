# Testing code

import numpy as np
import unittest
from Tectonic_Utils.geodesy import fault_vector_functions
from Tectonic_Utils.geodesy import insar_vector_functions

class Tests(unittest.TestCase):

    def test_strike(self):
        strike = fault_vector_functions.get_strike(deltax=1.0, deltay=0.0);
        self.assertEqual(strike, 90);
        strike = fault_vector_functions.get_strike(deltax=-1.0, deltay=0.0);
        self.assertEqual(strike, 270);
        strike = fault_vector_functions.get_strike(deltax=0.0, deltay=-1.0);
        self.assertEqual(strike, 180);
        angle = -160;
        strike = fault_vector_functions.get_strike(deltax=np.cos(np.deg2rad(angle)), deltay=np.sin(np.deg2rad(angle)));
        self.assertEqual(strike, 250);
        return;

    def test_rake(self):
        rake = fault_vector_functions.get_rake(strike_slip=1.0, dip_slip=0.5);
        self.assertAlmostEqual(rake, 26.5650511770779);
        return;

    def test_plane_normal(self):
        plane_normal = fault_vector_functions.get_plane_normal(strike=0, dip=0);
        self.assertAlmostEqual(plane_normal[0], 0);
        self.assertAlmostEqual(plane_normal[1], 0);
        self.assertAlmostEqual(plane_normal[2], 1);
        plane_normal = fault_vector_functions.get_plane_normal(strike=90, dip=89.99);
        self.assertAlmostEqual(plane_normal[0], 0);
        self.assertAlmostEqual(plane_normal[1], -1);
        self.assertAlmostEqual(plane_normal[2], 0.0001745329);
        plane_normal = fault_vector_functions.get_plane_normal(strike=270, dip=89.99);
        self.assertAlmostEqual(plane_normal[0], 0);
        self.assertAlmostEqual(plane_normal[1], 1);
        self.assertAlmostEqual(plane_normal[2], 0.0001745329);
        plane_normal = fault_vector_functions.get_plane_normal(strike=180, dip=1);
        self.assertAlmostEqual(plane_normal[0],  -0.0174524064372835);
        self.assertAlmostEqual(plane_normal[1], 0);
        self.assertAlmostEqual(plane_normal[2], 0.999847695156391);
        return;

    def test_lkv_conversion(self):
        """ Testing that look vectors are converted correctly"""
        [lkv_e, lkv_n, lkv_u] = insar_vector_functions.flight_incidence_angles2look_vector(190, 30);
        [out_flight, out_inc] = insar_vector_functions.look_vector2flight_incidence_angles(lkv_e, lkv_n, lkv_u);
        self.assertAlmostEqual(190, out_flight);
        self.assertAlmostEqual(30, out_inc);
        return;


if __name__ == "__main__":
    unittest.main();
