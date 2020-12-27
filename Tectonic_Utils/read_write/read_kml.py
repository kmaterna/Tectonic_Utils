# The purpose of this code is to read a simple box drawn in Google Earth and saved as a KML file. 
# The return value is a list of latitude and longitude coordinates. 


def read_simple_kml(infile):
	print("Reading %s into arrays of lon and lat..." % infile)
	start = 0;
	lats, lons = [], [];
	ifile = open(infile, 'r');
	for line in ifile:
		if start == 1:
			temp = line.split()
			for item in temp:
				lons.append(float(item.split(',')[0]))
				lats.append(float(item.split(',')[1]))
			break;
		if "coordinates" in line:
			start = 1;
	return lons, lats;


if __name__ == "__main__":
	ifile = "test1.kml"
	lons, lats = read_simple_kml(ifile);
