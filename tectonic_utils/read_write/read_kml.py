
from fastkml import kml


def read_simple_kml(infile):
    """
    Read a simple box drawn in Google Earth and saved as a KML file with field 'coordinates'.
    For more complicated KMLs with many elements, I'm switching to use the Python library "fastkml" for reading.

    :param infile: kml file with simple box
    :type infile: string
    :returns: lons, lats as lists that represent the coordinates of box vertices
    :rtype: list, list
    """
    print("Reading %s into arrays of lon and lat..." % infile)
    start = 0
    lats, lons = [], []
    ifile = open(infile, 'r')
    for line in ifile:
        if start == 1:
            temp = line.split()
            for item in temp:
                lons.append(float(item.split(',')[0]))
                lats.append(float(item.split(',')[1]))
            break
        if "coordinates" in line:
            start = 1
    return lons, lats


def read_kml_geometries_into_coord_lists(kml_name):
    """
    Iteratively traverse the etree inside a kml file, looking for things that have 'geometry'.
    This includes Placemarks: i.e., Point, LineString, or Polygon objects
    If they do have 'geometry', then we can print their coordinates and return them as lists.
    An attempt to create a GMT-compatible fault format from a KML

    :param kml_name: string, name of file
    :return: list of tuples or lists, containing coordinates of the drawn points
    """
    coord_lists = []

    with open(kml_name, "rb") as f:
        doc = f.read()

    k = kml.KML()
    k.from_string(doc)

    def walk_feature(feature, depth=0):
        indent = '  ' * depth
        print(f"{indent}- {type(feature).__name__}: {getattr(feature, 'name', None)}")

        if hasattr(feature, 'geometry') and feature.geometry:
            # print(f"{indent}  Geometry: {feature.geometry}")
            geom = feature.geometry
            coords = list(geom.coords) if hasattr(geom, 'coords') else geom
            # print(f"{indent}  Coordinates: {coords}")
            coord_lists.append(coords)

        if hasattr(feature, 'features'):
            for subfeature in feature.features():
                walk_feature(subfeature, depth + 1)

    # Start from the root object
    for feature in k.features():
        walk_feature(feature)

    return coord_lists
