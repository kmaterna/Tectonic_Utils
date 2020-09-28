# The functions in this file are used for reading common file types into structures in Pythin
# Example: a multi-segment file with polygons or lines, as could be used in GMT



def read_gmt_multisegment_latlon(fields_file, split_delimiter=' '):
	# GENERALIZED GMT MULTISEGMENT FILE READER IN PYTHON
    # Returns a list of lists, each one with a single segment
    print("reading gmt multisegment file %s" % fields_file);
    ifile = open(fields_file);
    lon_collection = [];
    lat_collection = [];
    lon_temp = [];
    lat_temp = [];
    for line in ifile:
        if line.split()[0] == '>>' or line.split()[0] == '>':
            if lon_temp != []:
                lon_collection.append(lon_temp);
                lat_collection.append(lat_temp);
            lon_temp = [];
            lat_temp = [];
            continue;
        else:
            temp = line.split(split_delimiter);
            lon_temp.append(float(temp[0]));
            lat_temp.append(float(temp[1]));
    lon_collection.append(lon_temp);
    lat_collection.append(lat_temp);
    return lon_collection, lat_collection;

