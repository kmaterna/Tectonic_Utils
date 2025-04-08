"""
Functions to read common file types into structures in Python.
Example: a multi-segment file with polygons or lines, as could be used in GMT.
"""


def read_gmt_multisegment_latlon(input_file, split_delimiter=' '):
    """
    Generalized GMT multisegment file reader.
    Returns lon and lat in a list of lists, each element with a single segment.

    :param input_file: name of input file
    :type input_file: string
    :param split_delimiter: delimiter between values on the same line, defaults to space
    :type split_delimiter: string, optional
    :returns: list of lons, list of lats
    :rtype: list
    """
    print("reading gmt multisegment file %s" % input_file)
    ifile = open(input_file)
    lon_collection, lat_collection = [], []
    lon_temp, lat_temp = [], []
    for line in ifile:
        if line.split()[0] == '>>' or line.split()[0] == '>':
            if lon_temp:
                lon_collection.append(lon_temp)
                lat_collection.append(lat_temp)
            lon_temp, lat_temp = [], []
            continue
        else:
            temp = line.split(split_delimiter)
            lon_temp.append(float(temp[0]))
            lat_temp.append(float(temp[1]))
    lon_collection.append(lon_temp)
    lat_collection.append(lat_temp)
    return lon_collection, lat_collection


def write_gmt_multisegment(list_of_coords_segments, filename):
    """
    Write a list of lists of coordiantes into a GMT-compatible multi-segment file
     [(-115.5650522767964, 33.11974272741214, 0.0),
     (-115.5642209309202, 33.12270979703938, 0.0),
     (-115.5637985114591, 33.12561963960839, 0.0)]

    :param list_of_coords_segments: list of lists, or list of tuples
    :param filename: string, name of output file
    """
    print("Writing file %s " % filename)
    with open(filename, 'w') as ofile:
        for segment in list_of_coords_segments:
            ofile.write('>\n')
            for coordinate in segment:
                ofile.write("%f %f \n" % (coordinate[0], coordinate[1]))
    return
