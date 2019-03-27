# Add a vector of km to a set of coordinate points

import numpy as np 


def add_vector_to_coords(lon0, lat0, dx, dy):
	lat1=lat0+dy/111.000;
	lon1=lon0+dx/(111.000*np.cos(np.deg2rad(lat0)));
	return_coords = [lon1, lat1]
	return return_coords;


if __name__=="__main__":
	lon0 = -124.45;
	lat0 = 40.305;
	top_left = add_vector_to_coords(lon0, lat0, -14.5, 3.0);
	top_right = add_vector_to_coords(lon0, lat0, 15.0, 3.0);
	bottom_left = add_vector_to_coords(lon0, lat0, -14.5, -3.5);
	bottom_right = add_vector_to_coords(lon0, lat0, 15.0, -3.5);

	print(top_left);
	print(top_right);
	print(bottom_left);
	print(bottom_right);
