

def wrap_lon(longitude):
    """
    :param longitude: float
    :returns: wrapped longitude, from -180 to 180
    """
    if longitude > 360 or longitude < -360:
        raise ValueError("longitude outside normal range (-360 to 360)");
    if longitude > 180:
        longitude = longitude - 360;
    if longitude < -180:
        longitude = longitude + 360;
    return longitude;
