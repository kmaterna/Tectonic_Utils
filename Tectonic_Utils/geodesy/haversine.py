#!/usr/bin/env python

# Haversine formula example in Python
# Author: Wayne Dyck
# From here, you can also add a vector of km lengths to a lon/lat coordinate, returning coordinates in lon/lat. 


import numpy as np 
import math

def distance(origin, destination):
    """
    Computes the distance [in km] between origin [lat1, lon1] and destination [lat2, lon2]. 
    Lat-lon must be specified in that order. 
    """
    lat1, lon1 = origin
    lat2, lon2 = destination
    radius = 6371  # km

    dlat = math.radians(lat2-lat1)
    dlon = math.radians(lon2-lon1)
    a = math.sin(dlat/2) * math.sin(dlat/2) + math.cos(math.radians(lat1)) \
        * math.cos(math.radians(lat2)) * math.sin(dlon/2) * math.sin(dlon/2)
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1-a))
    d = radius * c

    return d


def calculate_initial_compass_bearing(pointA, pointB):
    """
    Calculates the bearing between two points.
    The formulae used is the following:
        theta = atan2(sin(delta_long).cos(lat2),cos(lat1).sin(lat2) - sin(lat1).cos(lat2).cos(delta_long))
    :Parameters:
      - `pointA: The tuple representing the latitude/longitude for the
        first point. Latitude and longitude must be in decimal degrees
      - `pointB: The tuple representing the latitude/longitude for the
        second point. Latitude and longitude must be in decimal degrees
    :Returns:
      The bearing in degrees (CW from north, just like strike)
    :Returns Type:
      float
    """
    if (type(pointA) != tuple) or (type(pointB) != tuple):
        raise TypeError("Only tuples are supported as arguments")

    lat1 = math.radians(pointA[0])
    lat2 = math.radians(pointB[0])

    diffLong = math.radians(pointB[1] - pointA[1])

    x = math.sin(diffLong) * math.cos(lat2)
    y = math.cos(lat1) * math.sin(lat2) - (math.sin(lat1) * math.cos(lat2) * math.cos(diffLong))

    initial_bearing = math.atan2(x, y)

    # Now we have the initial bearing but math.atan2 return values
    # from -180 to + 180 which is not what we want for a compass bearing
    # The solution is to normalize the initial bearing as shown below
    initial_bearing = math.degrees(initial_bearing)
    compass_bearing = (initial_bearing + 360) % 360

    return compass_bearing

def calculate_endpoint_given_bearing(origin, bearing, angular_distance_degrees):
    # head a certain angular distance (degrees) along a bearing
    # starting from a point on earth.  Where do you end up? 
    # Origin is a two-vector that contains [latitude, longitude]
    # bearing is clockwise from north in degrees, like strike
    # phi2, lambda2 are the latitude and longitude of the destination point. 
    # theta is bearing in radians. 
    # delta is the angular distance in radians, d/R (d = distance, R = radius of Earth)
    # Source: https://www.movable-type.co.uk/scripts/latlong.html
    lat0 = origin[0];
    lon0 = origin[1];
    phi1 = np.deg2rad(lat0)
    lambda1 = np.deg2rad(lon0);
    delta = np.deg2rad(angular_distance_degrees);
    theta = np.deg2rad(bearing);
    phi2 = np.arcsin((np.sin(phi1)*np.cos(delta)) + (np.cos(phi1)*np.sin(delta)*np.cos(theta)));
    lambda2 = lambda1 + np.arctan2(np.sin(theta)*np.sin(delta)*np.cos(phi1), np.cos(delta)-np.sin(phi1)*np.sin(phi2));

    lat2 = np.rad2deg(phi2)
    lon2 = np.rad2deg(lambda2);  
    destination = [lat2, lon2]
    return destination; 


def xy_distance(ref_loc, sta_loc):
    """
    Distance between two latitude/longitude pairs, 
    given in x-distance and y-distance in meters
    (assuming flat surface between the points)
    Returns x and y in meters. 
    """
    radius = distance(ref_loc, sta_loc);
    bearing = calculate_initial_compass_bearing((ref_loc[0], ref_loc[1]), (sta_loc[0], sta_loc[1]))
    azimuth = 90 - bearing;
    x = radius * np.cos(np.deg2rad(azimuth)) * 1000;
    y = radius * np.sin(np.deg2rad(azimuth)) * 1000;
    return [x, y];


def add_vector_to_coords(lon0, lat0, dx, dy):
    # Add a vector of km to a set of latitude/longitude points. 
    lat1 = lat0+dy/111.000;
    lon1 = lon0+dx/(111.000*np.cos(np.deg2rad(lat0)));
    return_coords = [lon1, lat1]
    return return_coords;
