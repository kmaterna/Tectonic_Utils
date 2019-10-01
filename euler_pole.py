"""
A function to rotate 
a point by a known euler pole. 
"""

import numpy as np
import sys

def point_rotation_by_Euler_Pole(Point, Euler_Pole):

	R_point = get_r(Point[0], Point[1]);
	R_ep = get_r(Euler_Pole[0], Euler_Pole[1]);
	unit_ep = get_unit_vector(R_ep);
	omega_raw = degma2radyr(Euler_Pole[2]);
	omega = omega_raw * unit_ep; # in radians per year

	velocity_of_transformation = np.cross(omega, R_point);  # velocity at the station from the euler pole rotation
	velocity_of_transformation = velocity_of_transformation*1000;  # mm/yr in x, y, z
			
	xvel = velocity_of_transformation[0];
	yvel = velocity_of_transformation[1];
	zvel = velocity_of_transformation[2];
	[east_transform, north_transform, up_transform] = xyz2enu(xvel, yvel, zvel, Point[0], Point[1]);

	return [east_transform, north_transform, up_transform];



def get_unit_vector(vec):
	mag=np.sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]);
	vec=np.divide(vec,mag);
	return vec;

def degma2radyr(omega):
	# degrees/Ma to radians/yr
	radyr = omega*(np.pi/180)*1e-6;
	return radyr;

def get_r(lon,lat):
	# Vector from center of earth to the point in question
	# Definitions: 
	# 0 longitude: x = 0
	# 0 latitude: z = 0
	# North is positive z. 
	R_fixed = 6378000;  # In meters
	R_equatorial_disk = R_fixed * np.cos(np.deg2rad(lat)); 
	T_equatorial_disk = np.deg2rad(lon); 
	X = R_equatorial_disk*np.cos(T_equatorial_disk); 
	Y = R_equatorial_disk*np.sin(T_equatorial_disk); 
	Z = np.sqrt(R_fixed*R_fixed - X*X - Y*Y); 
	if lat<0:
		Z=Z*-1;
	return [X, Y, Z];

def get_unit_east(lon, lat):
	T_equatorial_disk = np.deg2rad(lon); 
	x = -np.sin(T_equatorial_disk);
	y = np.cos(T_equatorial_disk);
	return [x, y, 0];

def xyz2enu(x, y, z, lon, lat):
	vel_vector = [x, y, z];
	# Convert velocities from xyz to east north up
	# Assuming spherical earth
	# Assuming horizontal movement only 
	# The unit east vector is the only one we need. 
	# Then we take the dot product with the unit east vector. 
	# The north component is the remainder. 
	unit_east = get_unit_east(lon, lat);
	e = np.dot(vel_vector,unit_east);
	n = np.sqrt(x*x + y*y + z*z - e*e);
	if z < 0:
		n=n*-1;
	u = 0;
	return [e, n, u];


if __name__=="__main__":
	Euler_Pole = [69.9, -12.3, 0.55]; # Lon, Lat, Deg/Ma
	Point = [-124, 40.5];  # Lon, Lat
	[east_transform, north_transform, up_transform] = point_rotation_by_Euler_Pole(Point, Euler_Pole);
	total = np.sqrt(east_transform*east_transform + north_transform*north_transform);
	print("%.2f east, %.2f north, %.2f up, %.2f mm/yr total" % (east_transform, north_transform, up_transform, total) );
	

