"""
Useful utilities for defining fault planes and coordinate systems
"""

import numpy as np
import math
from Tectonic_Utils.geodesy import haversine


def xy2lonlat_single(xi, yi, reflon, reflat):
    lat = reflat + (yi * 1 / 111.000);
    lon = reflon + (xi * 1 / (111.000 * abs(np.cos(np.deg2rad(reflat)))));
    return lon, lat;


def latlon2xy_single(loni, lati, lon0, lat0):
    """returns the distance between a point and a reference in km."""
    radius = haversine.distance([lat0, lon0], [lati, loni]);
    bearing = haversine.calculate_initial_compass_bearing((lat0, lon0), (lati, loni))
    azimuth = 90 - bearing;
    x = radius * np.cos(np.deg2rad(azimuth));
    y = radius * np.sin(np.deg2rad(azimuth));
    return x, y;


def xy2lonlat(xi, yi, reflon, reflat):
    """can take a list or a single value"""
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


def latlon2xy(loni, lati, lon0, lat0):
    """can take a list or a single value"""
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


def get_plane_normal(strike, dip):
    """
    Get the unit normal to a plane of specified strike and dip. We first find the orthogonal unit vectors aligned
    with strike and dip directions that sit within the plane. The plane normal is their cross product,
    i.e. the outward facing unit normal vector, dip-cross-strike, in x-y-z coordinates.

    :param strike: strike, degrees
    :type strike: float
    :param dip: dip, degrees
    :type dip: float
    :return: 3-component vector of outward facing unit normal vector
    :rtype: np.ndarray
    """
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
    """Returns a unit vector in x-y-z coordinates"""
    theta = np.deg2rad(90 - strike);
    strike_vector = [np.cos(theta), np.sin(theta), 0];
    return strike_vector;


def get_dip_vector(strike, dip):
    """Returns a unit vector in x-y-z coordinates"""
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
    """Returns the strike of a line (in cw degrees from north) given the deltax and deltay in km."""
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
    """Just the pythagorean theorem"""
    return np.sqrt(strike_slip * strike_slip + dip_slip * dip_slip);


def get_strike_length(x0, x1, y0, y1):
    """Just the pythagorean theorem"""
    length = np.sqrt((x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0));
    return length;


def get_downdip_width(top, bottom, dip):
    W = abs(top - bottom) / np.sin(np.deg2rad(dip));  # guaranteed to be between 0 and 90
    return W;


def get_top_bottom_from_center(center_depth, width, dip):
    """Given a fault, where is the top and bottom?
    Width is total downdip width of the fault."""
    top = center_depth - (width / 2.0 * np.sin(np.deg2rad(dip)));
    bottom = center_depth + (width / 2.0 * np.sin(np.deg2rad(dip)));
    return top, bottom;


def get_top_bottom_from_top(top_depth, width, dip):
    bottom = top_depth + (width * np.sin(np.deg2rad(dip)));
    return top_depth, bottom;


def add_vector_to_point(x0, y0, vector_mag, vector_heading):
    """Vector heading defined as strike- CW from north."""
    theta = np.deg2rad(90 - vector_heading);
    x1 = x0 + vector_mag * np.cos(theta);
    y1 = y0 + vector_mag * np.sin(theta);
    return x1, y1;


def get_rake(strike_slip, dip_slip):
    """
    Positive slip is right lateral, and reverse.
    Range is -180 to 180.
    Will return 0 if dipslip,strikeslip==0,0
    """
    rake = np.rad2deg(math.atan2(dip_slip, strike_slip));
    return rake;
