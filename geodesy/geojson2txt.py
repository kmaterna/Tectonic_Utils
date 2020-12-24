# May 14, 2020
# In this script, I need utility functions to convert between 
# two different InSAR data representations
# 1. geojson, produced by Kite after downsampling
# 2. text reprsentation, used for Trever Hines' inversion code. 

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.patches import Polygon
import collections
import json

Downsampled_pixel = collections.namedtuple("Downsampled_pixel",
                                           ["mean", "median", "std", "BL_corner", "TR_corner", "unitE", "unitN",
                                            "unitU"]);
# Median and Mean values are in meters
# Unit vectors are from ground to satellite.


def read_geojson(infile):
    # This function reads a geojson as created by Kite downsampling
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


def plot_downsampled_InSAR(pixel_list, plot_name, vmin=None, vmax=None):
    plt.figure(figsize=(16, 8), dpi=300);
    color_boundary_object = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax, clip=True);
    custom_cmap = cm.ScalarMappable(norm=color_boundary_object, cmap='rainbow');

    # Getting lat and lon information for plotting.
    xdim = [pixel.BL_corner[0] for pixel in pixel_list];
    ydim = [pixel.BL_corner[1] for pixel in pixel_list];

    # Plotting each pixel
    for pixel in pixel_list:
        x_total = [pixel.BL_corner[0], pixel.BL_corner[0], pixel.TR_corner[0], pixel.TR_corner[0]]
        y_total = [pixel.BL_corner[1], pixel.TR_corner[1], pixel.TR_corner[1], pixel.BL_corner[1]]

        fault_vertices = np.column_stack((x_total, y_total));
        patch_color = custom_cmap.to_rgba(pixel.mean * 1000);  # in mm
        mypolygon = Polygon(fault_vertices, color=patch_color, alpha=1.0);
        ax = plt.gca();
        ax.add_patch(mypolygon);

    # Plot formatting
    plt.title('Downsampled %d pixels ' % (len(pixel_list)));
    plt.ylim([np.min(ydim) - 0.05, np.max(ydim) + 0.05]);
    plt.xlim([np.min(xdim) - 0.05, np.max(xdim) + 0.05]);
    plt.gca().tick_params(axis='both', which='major', labelsize=18)
    custom_cmap.set_array(np.arange(vmin, vmax, 100));
    cb = plt.colorbar(custom_cmap);
    cb.set_label('mm', fontsize=22);
    for l in cb.ax.yaxis.get_ticklabels():
        l.set_size(18)
    plt.xlabel("Longitude", fontsize=18);
    plt.ylabel("Latitude", fontsize=18);
    plt.savefig(plot_name);
    plt.close();

    # New figure for standard deviation of pixels
    plt.figure(figsize=(16, 8), dpi=300);
    color_boundary_object = matplotlib.colors.Normalize(vmin=0.0, vmax=40, clip=True);
    custom_cmap = cm.ScalarMappable(norm=color_boundary_object, cmap='rainbow');

    # Getting lat and lon information for plotting.
    xdim = [pixel.BL_corner[0] for pixel in pixel_list];
    ydim = [pixel.BL_corner[1] for pixel in pixel_list];

    # Plotting each pixel
    for pixel in pixel_list:
        x_total = [pixel.BL_corner[0], pixel.BL_corner[0], pixel.TR_corner[0], pixel.TR_corner[0]]
        y_total = [pixel.BL_corner[1], pixel.TR_corner[1], pixel.TR_corner[1], pixel.BL_corner[1]]

        fault_vertices = np.column_stack((x_total, y_total));
        patch_color = custom_cmap.to_rgba(pixel.std * 1000);  # in mm
        mypolygon = Polygon(fault_vertices, color=patch_color, alpha=1.0);
        ax = plt.gca();
        ax.add_patch(mypolygon);

    # Plot formatting
    plt.title('Downsampled %d pixels ' % (len(pixel_list)));
    plt.ylim([np.min(ydim) - 0.05, np.max(ydim) + 0.05]);
    plt.xlim([np.min(xdim) - 0.05, np.max(xdim) + 0.05]);
    plt.gca().tick_params(axis='both', which='major', labelsize=18)
    custom_cmap.set_array(np.arange(vmin, vmax, 100));
    cb = plt.colorbar(custom_cmap);
    cb.set_label('mm', fontsize=22);
    for l in cb.ax.yaxis.get_ticklabels():
        l.set_size(18)
    plt.xlabel("Longitude", fontsize=18);
    plt.ylabel("Latitude", fontsize=18);
    len_of_ext = len(plot_name.split('.')[-1]) + 1  # moving three spaces away from the end to add the file descriptor
    plt.savefig(plot_name[0:-len_of_ext] + "_std" + plot_name[-len_of_ext:]);
    plt.close();
    return;


def pixels_to_txt(pixel_list, text_file, bbox=[-180, 180, -90, 90], std_min=0.001):
    # Write in the format needed by Trever's inversion code
    # bbox is the optional bounding box with format [W,E,S,N];
    # std_min is the minimum value of uncertainty (m)
    # Lon Lat Value unitE unitN unitU
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
