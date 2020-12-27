import numpy as np 
from Tectonic_Utils.geodesy import euler_pole

if __name__=="__main__":
	Euler_Pole = [69.9, -12.3, 0.55]; # Lon, Lat, Deg/Ma
	Point = [-124, 40.5];  # Lon, Lat
	[east_transform, north_transform, up_transform] = euler_pole.point_rotation_by_Euler_Pole(Point, Euler_Pole);
	print("%.2f east, %.2f north, %.2f up" % (east_transform, north_transform, up_transform) );