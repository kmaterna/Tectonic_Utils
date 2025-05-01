"""
Haversine formula example in Python. Original author Wayne Dyck.
"""

import numpy as np 
import math


def distance(origin, destination, radius=6371):
    """
    Computes the distance between origin [lat1, lon1] and destination [lat2, lon2].

    :param origin: Tuple representing (latitude, longitude) of first point, in decimal degrees
    :type origin: array_like
    :param destination: Tuple representing (latitude, longitude) of second point, in decimal degrees
    :type destination: array_like
    :param radius: float, radius of Earth in km, default 6371
    :return: distance, in km
    :rtype: float
    """
    lat1, lon1 = origin
    lat2, lon2 = destination

    dlat = math.radians(lat2-lat1)
    dlon = math.radians(lon2-lon1)
    a = math.sin(dlat/2) * math.sin(dlat/2) + math.cos(math.radians(lat1)) \
        * math.cos(math.radians(lat2)) * math.sin(dlon/2) * math.sin(dlon/2)
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1-a))
    d = radius * c

    return d


def distance_vectorized(origin, destination, radius=6371.0):
    """
    Calculate the great-circle distance between pairs of points using numpy vectorized operations.

    Parameters:
    - origin: ndarray of shape (N, 2) with columns [lat, lon] in degrees
    - destination: ndarray of shape (N, 2) with columns [lat, lon] in degrees
    - radius: radius of the Earth in kilometers (default: 6371)

    Returns:
    - distances: ndarray of shape (N,) with distances in kilometers
    """
    lat1 = np.radians(origin[:, 0])
    lon1 = np.radians(origin[:, 1])
    lat2 = np.radians(destination[:, 0])
    lon2 = np.radians(destination[:, 1])

    dlat = lat2 - lat1
    dlon = lon2 - lon1

    a = np.sin(dlat / 2.0) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2.0) ** 2
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))

    dist = radius * c
    return dist


def calculate_initial_compass_bearing(pointA, pointB):
    r"""
    Calculate the bearing between two points.
    By the formula

    .. math::   \theta = atan2(\sin(\Delta_{long})*\cos(lat_2), \cos(lat_1)*\sin(lat_2) -
                         \sin(lat_1)*\cos(lat_2)*\cos(\Delta_{long})).

    :param pointA: Tuple representing (latitude, longitude) of first point, in decimal degrees
    :type pointA: array_like
    :param pointB: Tuple representing (latitude, longitude) of second point, in decimal degrees
    :type pointB: array_like
    :return: bearing, in degrees CW from north
    :rtype: float
    """
    if not isinstance(pointA, tuple) or not isinstance(pointB, tuple):
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


def calculate_initial_compass_bearing_vectorized(pointA, pointB):
    """
    Calculate initial compass bearing from pointA to pointB using numpy vectorized operations.

    Parameters:
    - pointA: ndarray of shape (N, 2) with [lat, lon] in degrees
    - pointB: ndarray of shape (N, 2) with [lat, lon] in degrees

    Returns:
    - bearings: ndarray of shape (N,) with compass bearings in degrees [0, 360)
    """
    lat1 = np.radians(pointA[:, 0])
    lat2 = np.radians(pointB[:, 0])
    diff_long = np.radians(pointB[:, 1] - pointA[:, 1])

    x = np.sin(diff_long) * np.cos(lat2)
    y = np.cos(lat1) * np.sin(lat2) - (
            np.sin(lat1) * np.cos(lat2) * np.cos(diff_long)
    )

    initial_bearing = np.arctan2(x, y)
    initial_bearing = np.degrees(initial_bearing)
    compass_bearing = (initial_bearing + 360) % 360

    return compass_bearing


def calculate_endpoint_given_bearing(origin, bearing, angular_distance_degrees):
    """
    Travel a certain angular distance (degrees) along a bearing starting from a coordinate.
    Source: https://www.movable-type.co.uk/scripts/latlong.html

    :param origin: Tuple representing (latitude, longitude) of first point, in decimal degrees
    :type origin: array_like
    :param bearing: angle clockwise from north, in decimal degrees
    :type bearing: float
    :param angular_distance_degrees: angular distance, in decimal degrees
    :type angular_distance_degrees: float
    :return: [lat, lon], in degrees
    :rtype: [float, float]
    """
    # Internally, phi2, lambda2 are the latitude and longitude of the destination point.
    # theta is bearing in radians.
    # delta is the angular distance in radians, d/R (d = distance, R = radius of Earth)

    lat0 = origin[0]
    lon0 = origin[1]
    phi1 = np.deg2rad(lat0)
    lambda1 = np.deg2rad(lon0)
    delta = np.deg2rad(angular_distance_degrees)
    theta = np.deg2rad(bearing)
    phi2 = np.arcsin((np.sin(phi1)*np.cos(delta)) + (np.cos(phi1)*np.sin(delta)*np.cos(theta)))
    lambda2 = lambda1 + np.arctan2(np.sin(theta)*np.sin(delta)*np.cos(phi1), np.cos(delta)-np.sin(phi1)*np.sin(phi2))
    lat2 = np.rad2deg(phi2)
    lon2 = np.rad2deg(lambda2)  
    destination = [lat2, lon2]
    return destination 


def xy_distance(ref_loc, sta_loc):
    """
    Pythagorean distance between two latitude/longitude pairs, assuming flat surface between the points.
    Returns x and y in meters.

    :param ref_loc: Tuple representing (latitude, longitude) of first point, in decimal degrees
    :type ref_loc: array_like
    :param sta_loc: Tuple representing (latitude, longitude) of second point, in decimal degrees
    :type sta_loc: array_like
    :return: [distance_x, distance_y], in m
    :rtype: [float, float]
    """
    radius = distance(ref_loc, sta_loc)  # in km
    bearing = calculate_initial_compass_bearing((ref_loc[0], ref_loc[1]), (sta_loc[0], sta_loc[1]))
    azimuth = 90 - bearing
    x = radius * np.cos(np.deg2rad(azimuth)) * 1000  # in m
    y = radius * np.sin(np.deg2rad(azimuth)) * 1000  # in m
    return [x, y]


def add_vector_to_coords(lon0, lat0, dx, dy):
    """Add a vector of km to a set of latitude/longitude points.

    :param lon0: Longitude of initial point, in degrees
    :type lon0: float
    :param lat0: Latitude of initial point, in degrees
    :type lat0: float
    :param dx: x component of vector added to the point, in km
    :type dx: float
    :param dy: y component of vector added to the point, in km
    :type dy: float
    :return: [lon1, lat1], in degrees
    :rtype: [float, float]
    """
    lat1 = lat0+dy/111.000
    lon1 = lon0+dx/(111.000*np.cos(np.deg2rad(lat0)))
    return [lon1, lat1]
