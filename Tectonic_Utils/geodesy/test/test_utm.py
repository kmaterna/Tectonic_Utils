import numpy as np
import unittest
from .. import utm_conversion


class Tests(unittest.TestCase):

    def test_utm(self):
        print("Testing utm conversion math.");
        lat = [40.3154333, 46.283900, 37.577833, 28.645650, 38.855550, 25.061783]
        lon = [-3.4857166, 7.8012333, -119.95525, -17.759533, -94.7990166, 121.640266]
        x, y, utm = utm_conversion.deg2utm(lat, lon);
        x_predicted = [458731.0, 407653.0, 239027.0, 230253.0, 343898.0, 362850.0];
        y_predicted = [4462881.0, 5126290.0, 4163083.0, 3171843.0, 4302285.0, 2772478.0];
        utm_predicted = ['30.0 T', '32.0 T', '11.0 S', '28.0 R', '15.0 S', '51.0 R'];
        self.assertListEqual(x_predicted, list(np.round(x)), 'deg2utm failed.');
        self.assertListEqual(y_predicted, list(np.round(y)), 'deg2utm failed.');
        self.assertListEqual(utm_predicted, list(utm), 'deg2utm failed.');

        x = [458731,  407653,  239027,  230253,  343898,  362850];
        y = [4462881, 5126290, 4163083, 3171843, 4302285, 2772478];
        utmzone = ['30 T', '32 T', '11 S', '28 R', '15 S', '51 R'];
        lat, lon = utm_conversion.utm2deg(x, y, utmzone);
        lat_predicted = [40.31543016, 46.28390185, 37.57783421, 28.64564715, 38.85555224, 25.06177977];
        lon_predicted = [-3.48571304, 7.80123464, -119.9552465, -17.75953665, -94.79901928, 121.64026588];
        self.assertListEqual(list(np.round(lat, 8)), lat_predicted, 'utm2deg failed.');
        self.assertListEqual(list(np.round(lon, 8)), lon_predicted, 'utm2deg failed.');
        return;


if __name__ == '__main__':
    unittest.main()
