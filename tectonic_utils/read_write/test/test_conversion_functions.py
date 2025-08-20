# Testing code

import numpy as np
import unittest
import os
from .. import netcdf_read_write as netrw


class Tests(unittest.TestCase):

    def test_pixel_node_writer(self):
        """
        See if the writing function for pixel-node files produces a pixel-node file.
        The behavior has been finicky for float32 vs float64
        Writing a full test for float32 would be good (although the example grd file gets pretty close)
        """
        grid_def = [-120, -114, 32, 37]
        inc = [0.02, 0.02]
        filename = 'test_outfile.nc'
        lons = np.arange(grid_def[0], grid_def[1] + 0.00001, inc[0])
        lats = np.arange(grid_def[2], grid_def[3] + 0.00001, inc[1])

        # Test a write function
        grid = np.zeros((len(lats), len(lons)))
        netrw.write_netcdf4(lons, lats, grid, filename)
        netrw.parse_pixelnode_registration(filename)
        os.remove(filename)
        os.remove('gmt.history')

        # Test a read-write cycle on an example grid
        [x, y, z] = netrw.read_any_grd(os.path.join("tectonic_utils", "read_write", "test", "example_grd.grd"))
        netrw.write_netcdf4(x, y, z, os.path.join("tectonic_utils", "read_write", "test", "written_example.grd"))
        netrw.parse_pixelnode_registration(os.path.join("tectonic_utils", "read_write", "test", "written_example.grd"))
        os.remove('gmt.history')

        return


if __name__ == "__main__":
    unittest.main()
