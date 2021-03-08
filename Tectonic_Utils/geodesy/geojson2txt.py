"""
This script contains utility functions to convert between two different InSAR data representations:

* geojson, produced by Kite after downsampling
* text representation, used for inversions
"""

import numpy as np
import collections
import json

Downsampled_pixel = collections.namedtuple("Downsampled_pixel",
                                           ["mean", "median", "std", "BL_corner", "TR_corner", "unitE", "unitN",
                                            "unitU"]);
Downsampled_pixel.__doc__ = """
:param mean: mean LOS value for pixel (meters)
:param median: median LOS value for pixel (meters)
:param std: standard deviation LOS value for pixel (meters)
:param BL_corner: Bottom Left corner coordinate (longitude, latitude)
:param TR_corner: Top Right corner coordinate (longitude, latitude)
:param unitE: east component of unit vector from ground to satellite
:param unitN: north component of unit vector from ground to satellite
:param unitU: up component of unit vector from ground to satellite
"""


def read_geojson(infile):
    """
    Reads a geojson as created by Kite downsampling into downsampled pixel objects.

    :param infile: name of geojson file
    :type infile: string
    :return: list of pixel objects
    :rtype: list
    """
    # "features" is a list with n elements
    # Each element is a downsampled pixel and is stored as a dictionary
    # containing 'type','id','geometry','properties'
    # These properties will be unpacked into a named tuple that contains helpful information.

    with open(infile) as f:
        data = json.load(f);
    features = data["features"];  # beginning to unpack the geojson
    pixel_list = [];
    for feature in features:
        # Unpacking what's in each pixel
        bl = feature["geometry"]["coordinates"][0][0];
        tr = feature["geometry"]["coordinates"][0][2];
        mean = feature["properties"]["mean"] * 0.001
        median = feature["properties"]["median"] * 0.001
        std = feature["properties"]["std"] * 0.001
        unitE = feature["properties"]["unitE"]
        unitN = feature["properties"]["unitN"]
        unitU = feature["properties"]["unitU"]
        onePixel = Downsampled_pixel(mean=mean, median=median, std=std, BL_corner=bl, TR_corner=tr,
                                     unitE=unitE, unitN=unitN, unitU=unitU);
        pixel_list.append(onePixel);
    return pixel_list;


def pixels_to_txt(pixel_list, text_file, bbox=(-180, 180, -90, 90), std_min=0.001):
    """
    Writes InSAR pixels in the format needed by Slippy inversion code: Lon Lat Value unitE unitN unitU

    :param pixel_list: list of pixel objects
    :type pixel_list: list
    :param text_file: name of output file
    :type text_file: str
    :param bbox: tuple of (W, E, S, N) bounding box, defaults to (-180, 180, -90, 90)
    :type bbox: array_like, optional
    :param std_min: minimum value of uncertainty (m), defaults to 0.001
    :type std_min: float, optional
    """
    ofile = open(text_file, 'w');
    ofile.write("# Header: lon, lat, disp(m), sig(m), unitE, unitN, unitU from ground to satellite\n")
    for pixel in pixel_list:
        lon = np.mean([pixel.BL_corner[0], pixel.TR_corner[0]]);
        lat = np.mean([pixel.BL_corner[1], pixel.TR_corner[1]]);
        std = np.max([pixel.std, std_min]);  # don't let uncertainty get unreasonably small
        if bbox[0] <= lon <= bbox[1]:
            if bbox[2] <= lat <= bbox[3]:
                ofile.write("%.5f %.5f %.5f %.5f %.5f %.5f %.5f\n" % (
                    lon, lat, pixel.median, std, pixel.unitE, pixel.unitN, pixel.unitU));
    ofile.close();
    print("Writing %s with %d pixels " % (text_file, len(pixel_list)));
    return;
