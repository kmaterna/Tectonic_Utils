"""
The functions in this script convert between fault formats
# Format 1: json format for Slippy
# Format 2: slip distribution format for Slippy
# Format 3: .intxt, Kathryn Materna's format for Elastic_stresses_py

The internal format here is a dictionary containing:
Format internal: strike(deg), dip(deg), length(km), width(km), lon(corner), lat(corner), depth(km), rake(deg), slip(m)
If the fault is a receiver fault, we put slip = 0

Format json: basis1, basis2, length(m), width(m), nlength, nwidth, strike, dip, position [lon, lat, dep], penalty

Format slippy: lon lat depth[m] strike[deg] dip[deg] length[m] width[m] left-lateral[m] thrust[m] tensile[m]

Format intxt: strike rake dip length(km) width(km) updip_corner_lon updip_corner_lat updip_corner_dep(km) slip
"""

import numpy as np
import json
from Tectonic_Utils.geodesy import haversine
import math


# ----------- IO FUNCTIONS ---------------- # 

def read_faults_intxt(infile):
    # Read all faults that are sources or receivers
    # Reads faults into a list of fault dictionaries
    # The lat/lon refer to the top left corner of the fault.
    fault_list = [];
    ifile = open(infile, 'r');
    for line in ifile:
        temp = line.split();
        if temp[0] == "S:" or temp[0] == "R:":
            one_fault = {"strike": float(temp[1]), "dip": float(temp[3]), "length": float(temp[4]),
                         "width": float(temp[5]), "lon": float(temp[6]), "lat": float(temp[7]), "depth": float(temp[8]),
                         "rake": float(temp[2])};
            if temp[0] == "R:":
                one_fault["slip"] = 0;
            else:
                one_fault["slip"] = float(temp[9]);
            fault_list.append(one_fault);
    ifile.close();
    return fault_list;


def read_faults_json(infile):
    # Read all faults from a json file (just geometry; don't have slip or rake)
    # It has to convert from fault center to fault corner.
    # Reads faults into a list of fault dictionaries
    # Faults read from JSON have zero slip
    fault_list = [];
    config_file = open(infile, 'r')
    config = json.load(config_file);
    for key in config["faults"].keys():
        one_fault = {"strike": config["faults"][key]["strike"], "dip": config["faults"][key]["dip"],
                     "length": config["faults"][key]["length"] / 1000.0,
                     "width": config["faults"][key]["width"] / 1000.0};
        center_lon = config["faults"][key]["position"][0];
        center_lat = config["faults"][key]["position"][1];
        x_start, y_start = add_vector_to_point(0, 0, one_fault["length"] / 2, one_fault["strike"] - 180);  # in km
        corner_lon, corner_lat = xy2lonlat(x_start, y_start, center_lon, center_lat);  #
        one_fault["lon"] = corner_lon;
        one_fault["lat"] = corner_lat;
        one_fault["depth"] = -config["faults"][key]["position"][2] / 1000;
        one_fault["rake"] = 0;
        one_fault["slip"] = 0;
        fault_list.append(one_fault);
    return fault_list;


def read_slippy_distribution(infile):
    # Read a file from the Slippy inversion outputs lon[degrees] lat[degrees] depth[m] strike[degrees] dip[degrees]
    # length[m] width[m] left-lateral[m] thrust[m] tensile[m] segment_num.
    # Lon/lat usually refer to the center top of the fault.
    # Must convert the lon/lat to the top left corner.
    fault_list = [];
    [lon, lat, depth, strike, dip, length, width, ll_slip,
     thrust_slip, _, _] = np.loadtxt(infile, skiprows=1, unpack=True, dtype={"names": ('lon', 'lat', 'depth', 'strike',
                                                                                       'dip', 'length', 'width', 'ss',
                                                                                       'ds', 'tensile', 'num'),
                                                                             "formats": (np.float, np.float, np.float,
                                                                                         np.float, np.float, np.float,
                                                                                         np.float, np.float, np.float,
                                                                                         np.float, np.float)});
    for i in range(len(lon)):
        one_fault = {"strike": strike[i], "dip": dip[i], "length": length[i] / 1000, "width": width[i] / 1000,
                     "depth": -depth[i] / 1000};
        center_lon = lon[i];
        center_lat = lat[i];
        x_start, y_start = add_vector_to_point(0, 0, one_fault["length"] / 2, one_fault["strike"] - 180);  # in km
        corner_lon, corner_lat = xy2lonlat(x_start, y_start, center_lon, center_lat);
        one_fault["lon"] = corner_lon;
        one_fault["lat"] = corner_lat;
        one_fault["rake"] = get_rake(ll_slip[i], thrust_slip[i]);
        one_fault["slip"] = get_total_slip(ll_slip[i], thrust_slip[i]);
        fault_list.append(one_fault);
    return fault_list;


def write_faults_intxt(faults, outfile, receiver=True, source=False, write_header=False):
    # Writes the files as receivers with zero slip, or sources with finite slip.
    ofile = open(outfile, 'w');
    ofile.write(
        '#S/R: strike rake dip length(km) width(km) updip_corner_lon updip_corner_lat updip_corner_dep(km) slip\n')
    if write_header:
        ofile.write("G: 0.250 0.40 -115.75 -115.25 -115.5 32.75 33.25 33.0 # general information\n");
        ofile.write("R: 191 -103 45 6 7 -115.507 33.1 0  # random fault\n");
    for fault in faults:
        if receiver is True:
            ofile.write("R: %.1f nan %.1f %.2f %.2f %.3f %.3f %.2f \n" % (
                fault["strike"], fault["dip"], fault["length"], fault["width"], fault["lon"], fault["lat"],
                fault["depth"]));
        if source is True:
            ofile.write("S: %.1f %.1f %.1f %.2f %.2f %.3f %.3f %.2f %.3f \n" % (
                fault["strike"], fault["rake"], fault["dip"], fault["length"], fault["width"], fault["lon"],
                fault["lat"],
                fault["depth"], fault["slip"]));
    ofile.close();
    return;


def write_faults_json(faults_list, outfile):
    output = {};
    faults = {};
    count = 0;
    for fault in faults_list:
        # Convert the fault (which has top left corner) into a fault with top center coordinate
        corner_lon = fault["lon"]
        corner_lat = fault["lat"]
        x_center, y_center = add_vector_to_point(0, 0, fault["length"] / 2, fault["strike"]);  # in km
        center_lon, center_lat = xy2lonlat(x_center, y_center, corner_lon, corner_lat);  #

        count = count + 1;
        label = "fault" + str(count)
        fault["length"] = fault["length"] * 1000;
        fault["width"] = fault["width"] * 1000;
        fault["basis1"] = [1, 0, 0];
        fault["basis2"] = None;
        fault["Nlength"] = 1;
        fault["Nwidth"] = 1;
        fault["penalty"] = 1;
        fault["position"] = [center_lon, center_lat, fault["depth"]];
        fault.pop("lon");
        fault.pop("lat");
        fault.pop("depth");
        faults[label] = fault;
    output["faults"] = faults;
    output["plotter"] = "gmt";
    with open(outfile, 'w') as ofile:
        json.dump(output, ofile, indent=4);
    return;


def write_slipdistribution(faults, outfile):
    # Write a slip distribution for Slippy.
    # lon[deg] lat[deg] depth[m] strike[deg] dip[deg] len[m] wid[m] left-lateral[m] thrust[m] tensile[m] segment_num
    # Remember depths are negative in meters
    # Remeber lon/lat refer to fault center.
    # Will write this function later.
    ofile = open(outfile, 'w');
    ofile.write(
        '# lon[degrees] lat[degrees] depth[m] strike[degrees] dip[degrees] length[m] width[m] left-lateral[m] thrust['
        'm] tensile[m] segment_num\n');
    # for fault in faults:
    print("THIS FUNCTION ISN'T WRITTEN YET");
    ofile.close();
    return;


# ------- CONVERSION FUNCTIONS ------------ # 

def json2intxt(infile, outfile):
    faults = read_faults_json(infile);
    write_faults_intxt(faults, outfile, receiver=True, source=False);
    return;


def intxt2json(infile, outfile):
    faults = read_faults_intxt(infile);
    write_faults_json(faults, outfile);
    return;


def slippydist2intxt(infile, outfile):
    faults = read_slippy_distribution(infile);
    write_faults_intxt(faults, outfile, receiver=False, source=True, write_header=True);
    return;


def intxt2slipdistribution(infile, outfile):
    faults = read_faults_intxt(infile);
    write_slipdistribution(faults, outfile);
    return;


# ------------- MATH FUNCTIONS ------------------ # 


def get_plane_normal(strike, dip):
    # Given a strike and dip, find orthogonal unit vectors
    # aligned with strike and dip directions that sit within the plane.
    # The plane normal is their cross product.
    # Returns in x, y, z coordinates.
    strike_vector = get_strike_vector(strike);  # unit vector
    dip_vector = get_dip_vector(strike, dip);  # unit vector
    plane_normal = np.cross(dip_vector, strike_vector);  # dip x strike for outward facing normal, by right hand rule.
    return plane_normal;


def get_dip_degrees(x0, y0, z0, x1, y1, z1):
    horizontal_length = get_strike_length(x0, x1, y0, y1);
    vertical_distance = abs(z1 - z0);
    dip = np.rad2deg(math.atan2(vertical_distance, horizontal_length));
    return dip;


def get_strike_vector(strike):
    # Returns a unit vector in x-y-z coordinates
    theta = np.deg2rad(90 - strike);
    strike_vector = [np.cos(theta), np.sin(theta), 0];
    return strike_vector;


def get_dip_vector(strike, dip):
    # Returns a unit vector in x-y-z coordinates
    downdip_direction_theta = np.deg2rad(-strike);  # theta(strike+90)
    dip_unit_vector_z = np.sin(np.deg2rad(dip))  # the vertical component of the downdip unit vector
    dip_unit_vector_xy = np.sqrt(
        1 - dip_unit_vector_z * dip_unit_vector_z);  # the horizontal component of the downdip unit vector
    dip_vector = [dip_unit_vector_xy * np.cos(downdip_direction_theta),
                  dip_unit_vector_xy * np.sin(downdip_direction_theta), -dip_unit_vector_z];
    return dip_vector;


def get_vector_magnitude(vector):
    magnitude = 0;
    total = 0;
    for i in range(len(vector)):
        total = total + vector[i] * vector[i];
        magnitude = np.sqrt(total);
    return magnitude;


def get_strike(deltax, deltay):
    # Returns the strike of a line (in cw degrees from north) given the deltax and deltay in km.
    slope = math.atan2(deltay, deltax);
    strike = 90 - np.rad2deg(slope);
    if strike < 0:
        strike = strike + 360;
    return strike;


def get_rtlat_dip_slip(slip, rake):
    strike_slip = -slip * np.cos(np.deg2rad(rake));  # negative sign for convention of right lateral slip
    dip_slip = slip * np.sin(np.deg2rad(rake));
    return strike_slip, dip_slip;


def get_total_slip(strike_slip, dip_slip):
    # Just the pythagorean theorem
    return np.sqrt(strike_slip * strike_slip + dip_slip * dip_slip);


def get_strike_length(x0, x1, y0, y1):
    # Just the pythagorean theorem
    length = np.sqrt((x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0));
    return length;


def get_downdip_width(top, bottom, dip):
    W = abs(top - bottom) / np.sin(np.deg2rad(dip));  # guaranteed to be between 0 and 90
    return W;


def get_top_bottom(center_depth, width, dip):
    # Given a fault, where is the top and bottom?
    # Width is total downdip width of the fault.
    top = center_depth - (width / 2.0 * np.sin(np.deg2rad(dip)));
    bottom = center_depth + (width / 2.0 * np.sin(np.deg2rad(dip)));
    return top, bottom;


def get_top_bottom_from_top(top_depth, width, dip):
    bottom = top_depth + (width * np.sin(np.deg2rad(dip)));
    return top_depth, bottom;


def add_vector_to_point(x0, y0, vector_mag, vector_heading):
    # Vector heading defined as strike- CW from north.
    theta = np.deg2rad(90 - vector_heading);
    x1 = x0 + vector_mag * np.cos(theta);
    y1 = y0 + vector_mag * np.sin(theta);
    return x1, y1;


def get_rake(strike_slip, dip_slip):
    # Positive slip is right lateral, and reverse.
    # Range is -180 to 180.
    # Will return 0 if dipslip,strikeslip==0,0
    rake = np.rad2deg(math.atan2(dip_slip, strike_slip));
    return rake;


def get_fault_center(fault_object):
    # Compute the x-y-z coordinates of the center of a fault patch.
    # Index is the i'th fault patch in this fault_object
    W = get_downdip_width(fault_object.top, fault_object.bottom, fault_object.dipangle);
    center_z = (fault_object.top + fault_object.bottom) / 2.0;
    updip_center_x = (fault_object.xstart + fault_object.xfinish) / 2.0;
    updip_center_y = (fault_object.ystart + fault_object.yfinish) / 2.0;
    vector_mag = W * np.cos(
        np.deg2rad(fault_object.dipangle)) / 2.0;  # how far the middle is displaced downdip from map-view
    center_point = add_vector_to_point(updip_center_x, updip_center_y, vector_mag,
                                       fault_object.strike + 90);  # strike+90 = downdip direction.
    center = [center_point[0], center_point[1], center_z];
    return center;


def get_fault_four_corners(fault_object):
    # Get the four corners of the object, including updip and downdip.
    W = get_downdip_width(fault_object.top, fault_object.bottom, fault_object.dipangle);
    strike = fault_object.strike;

    updip_point0 = [fault_object.xstart, fault_object.ystart];
    updip_point1 = [fault_object.xfinish, fault_object.yfinish];
    vector_mag = W * np.cos(
        np.deg2rad(fault_object.dipangle));  # how far the bottom edge is displaced downdip from map-view
    downdip_point0 = add_vector_to_point(fault_object.xstart, fault_object.ystart, vector_mag,
                                         strike + 90);  # strike+90 = downdip direction.
    downdip_point1 = add_vector_to_point(fault_object.xfinish, fault_object.yfinish, vector_mag, strike + 90);

    x_total = [updip_point0[0], updip_point1[0], downdip_point1[0], downdip_point0[0], updip_point0[0]];
    y_total = [updip_point0[1], updip_point1[1], downdip_point1[1], downdip_point0[1], updip_point0[1]];
    x_updip = [updip_point0[0], updip_point1[0]];
    y_updip = [updip_point0[1], updip_point1[1]];
    return [x_total, y_total, x_updip, y_updip];


def xy2lonlat(xi, yi, reflon, reflat):  # can take a list or a single value
    if type(xi) == float or type(xi) == np.float64 or type(
            xi) == int:  # if we are getting a single value, we return a single value.
        lon, lat = xy2lonlat_single(xi, yi, reflon, reflat);
    else:  # if we are getting a list of values, we return a list of the same dimensions
        lat, lon = [], [];
        for i in range(len(xi)):
            loni, lati = xy2lonlat_single(xi[i], yi[i], reflon, reflat);
            lon.append(loni);
            lat.append(lati);
    return lon, lat;


def latlon2xy(loni, lati, lon0, lat0):  # can take a list or a single value
    if type(loni) == float or type(loni) == np.float64 or type(
            loni) == int:  # if we are getting a single value, we return a single value.
        x, y = latlon2xy_single(loni, lati, lon0, lat0);
    else:  # If we are getting a list, return a list of the same dimensions
        x, y = [], [];
        for i in range(len(loni)):
            xi, yi = latlon2xy_single(loni[i], lati[i], lon0, lat0);
            x.append(xi);
            y.append(yi);
    return [x, y];


# THE MATH FUNCTIONS
def xy2lonlat_single(xi, yi, reflon, reflat):
    lat = reflat + (yi * 1 / 111.000);
    lon = reflon + (xi * 1 / (111.000 * abs(np.cos(np.deg2rad(reflat)))));
    return lon, lat;


def latlon2xy_single(loni, lati, lon0, lat0):
    # returns the distance between a point and a reference in km.
    radius = haversine.distance([lat0, lon0], [lati, loni]);
    bearing = haversine.calculate_initial_compass_bearing((lat0, lon0), (lati, loni))
    azimuth = 90 - bearing;
    x = radius * np.cos(np.deg2rad(azimuth));
    y = radius * np.sin(np.deg2rad(azimuth));
    return x, y;
