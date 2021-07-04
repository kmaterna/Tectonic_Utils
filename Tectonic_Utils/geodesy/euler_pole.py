"""
Functions to rotate a point by a known euler pole.
"""

import numpy as np
from . import fault_vector_functions


def point_rotation_by_Euler_Pole(Point, Euler_Pole):
    """
    Compute the velocity of rotation of a point about an Euler pole on a spherical earth.
    This function is useful for computing the velocity of a stationary point in one reference frame
    with respect to another reference frame.
    The resulting velocity is assumed to be horizontal.

    :param Point: [longitude, latitude] of observation point, in degrees
    :type Point: array_like
    :param Euler_Pole: [longitude, latitude, omega] of Euler Pole, in degrees and degrees/Ma
    :type Euler_Pole: array_like
    :returns: [e_velocity, n_velocity, u_velocity] of point in rotated reference frame, in mm/yr
    :rtype: array_like
    """
    R_point = get_r(Point[0], Point[1]);
    R_ep = get_r(Euler_Pole[0], Euler_Pole[1]);
    unit_ep = fault_vector_functions.get_unit_vector(R_ep);
    omega_raw = degma2radyr(Euler_Pole[2]);
    omega = omega_raw * unit_ep;  # in radians per year

    velocity_of_transformation = np.cross(omega, R_point);  # velocity at the station from the euler pole rotation
    velocity_of_transformation = velocity_of_transformation * 1000;  # mm/yr in x, y, z

    xvel = velocity_of_transformation[0];
    yvel = velocity_of_transformation[1];
    zvel = velocity_of_transformation[2];
    [east_transform, north_transform] = xyz2en(xvel, yvel, zvel, Point[0]);
    up_transform = 0;  # by definition the velocity will be horizontal

    return [east_transform, north_transform, up_transform];


def degma2radyr(omega):
    """Convert omega from degrees/Ma to radians/yr"""
    radyr = omega * (np.pi / 180) * 1e-6;
    return radyr;


def get_r(lon, lat):
    """
    Vector from center of earth to the point in question assuming a spherical earth.
    The XYZ coordinate system has x=0 at longitude=0 and z=0 at the equator with positive to the north.

    :param lon: Longitude of initial point, in degrees
    :type lon: float
    :param lat: Latitude of initial point, in degrees
    :type lat: float
    :returns: [X, Y, Z] coordinates in meters.
    :rtype: [float, float, float]
    """
    R_fixed = 6378000;  # In meters
    R_equatorial_disk = R_fixed * np.cos(np.deg2rad(lat));
    T_equatorial_disk = np.deg2rad(lon);
    X = R_equatorial_disk * np.cos(T_equatorial_disk);
    Y = R_equatorial_disk * np.sin(T_equatorial_disk);
    Z = np.sqrt(R_fixed * R_fixed - X * X - Y * Y);
    if lat < 0:
        Z = Z * -1;
    return [X, Y, Z];


def get_unit_east(lon):
    """
    Unit east vector from a point on earth's surface in XYZ coordinates.
    The XYZ coordinate system has x=0 at longitude=0 and z=0 at the equator with positive to the north.
    The return value of Z is zero for eastward motion.

    :param lon: Longitude of initial point, in degrees
    :type lon: float
    :returns: [X, Y, Z] components
    :rtype: [float, float, float]
    """
    T_equatorial_disk = np.deg2rad(lon);
    x = -np.sin(T_equatorial_disk);
    y = np.cos(T_equatorial_disk);
    return [x, y, 0];


def xyz2en(x, y, z, lon):
    """
    Convert velocities from xyz to horizontal east and north, assuming spherical earth and no vertical motion.
    We take the dot product of the velocity with the unit east vector and the north component is the remainder.
    A more complex function xyz2enu(X, Y, Z, lon, lat) could be written later.

    :param x: x velocity at observation point
    :type x: float
    :param y: y velocity at observation point
    :type y: float
    :param z: z velocity at observation point
    :type z: float
    :param lon: Longitude of observation point, in degrees
    :type lon: float
    :returns: [east_vel, north_vel]
    :rtype: [float, float]
    """
    vel_vector = [x, y, z];
    unit_east = get_unit_east(lon);
    e = np.dot(vel_vector, unit_east);
    n = np.sqrt(x * x + y * y + z * z - e * e);
    if z < 0:
        n = n * -1;
    return [e, n];


if __name__ == "__main__":
    Euler_Pole = [69.9, -12.3, 0.55];  # Lon, Lat, Deg/Ma
    Point = [-124, 40.5];  # Lon, Lat
    [east_transform, north_transform, up_transform] = point_rotation_by_Euler_Pole(Point, Euler_Pole);
    total = np.sqrt(east_transform * east_transform + north_transform * north_transform);
    print("%.2f east, %.2f north, %.2f up, %.2f mm/yr total" % (east_transform, north_transform, up_transform, total));
