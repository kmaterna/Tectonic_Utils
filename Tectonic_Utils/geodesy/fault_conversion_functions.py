"""
The functions in this script convert between fault formats

* Format 1: json format for Slippy
* Format 2: slip distribution format for Slippy
* Format 3: .intxt, Kathryn Materna's format for Elastic_stresses_py

The internal format here is a dictionary containing:
Format internal: {strike(deg), dip(deg), length(km), width(km), lon(corner), lat(corner), depth(km), rake(deg), slip(m)}
If the fault is a receiver fault, we put slip = 0

Format json: basis1, basis2, length(m), width(m), nlength, nwidth, strike, dip, position [lon, lat, dep], penalty

Format slippy: lon lat depth[m] strike[deg] dip[deg] length[m] width[m] left-lateral[m] thrust[m] tensile[m]

Format intxt: strike rake dip length(km) width(km) updip_corner_lon updip_corner_lat updip_corner_dep(km) slip
"""

import numpy as np
import json
from Tectonic_Utils.geodesy import fault_vector_functions


def read_faults_intxt(infile):
    """
    Read all faults that are sources or receivers into a list of fault dictionaries.
    The lat/lon refer to the top left corner of the fault.

    :param infile: name of input .intxt file
    :type infile: string
    :returns: list of fault dictionaries
    :rtype: list
    """
    fault_list = [];
    ifile = open(infile, 'r');
    for line in ifile:
        temp = line.split();
        if temp[0] == "Source_Patch:" or temp[0] == "Receiver:":
            one_fault = {"strike": float(temp[1]), "dip": float(temp[3]), "length": float(temp[4]),
                         "width": float(temp[5]), "lon": float(temp[6]), "lat": float(temp[7]), "depth": float(temp[8]),
                         "rake": float(temp[2])};
            if temp[0] == "Receiver:":
                one_fault["slip"] = 0;
            else:
                one_fault["slip"] = float(temp[9]);
            fault_list.append(one_fault);
    ifile.close();
    return fault_list;


def read_faults_json(infile):
    """
    Read all faults from a json file (just geometry; no slip or rake) into a list of fault dictionaries.
    It has to convert from fault center to fault corner.
    Faults read from JSON have zero slip.

    :param infile: name of input json file
    :type infile: string
    :returns: list of fault dictionaries
    :rtype: list
    """
    fault_list = [];
    config_file = open(infile, 'r')
    config = json.load(config_file);
    for key in config["faults"].keys():
        one_fault = {"strike": config["faults"][key]["strike"], "dip": config["faults"][key]["dip"],
                     "length": config["faults"][key]["length"] / 1000.0,
                     "width": config["faults"][key]["width"] / 1000.0};
        center_lon = config["faults"][key]["position"][0];
        center_lat = config["faults"][key]["position"][1];
        x_start, y_start = fault_vector_functions.add_vector_to_point(0, 0, one_fault["length"] / 2,
                                                                      one_fault["strike"] - 180);  # in km
        corner_lon, corner_lat = fault_vector_functions.xy2lonlat(x_start, y_start, center_lon, center_lat);  #
        one_fault["lon"] = corner_lon;
        one_fault["lat"] = corner_lat;
        one_fault["depth"] = -config["faults"][key]["position"][2] / 1000;
        one_fault["rake"] = 0;
        one_fault["slip"] = 0;
        fault_list.append(one_fault);
    return fault_list;


def read_slippy_distribution(infile):
    """
    Read a file from the Slippy inversion outputs (lon[degrees] lat[degrees] depth[m] strike[degrees] dip[degrees]
    length[m] width[m] left-lateral[m] thrust[m] tensile[m] segment_num) into a list of fault dictionaries.
    Lon/lat usually refer to the center top of the fault, so it must convert the lon/lat to the top left corner.

    :param infile: name of input slip distribution file
    :type infile: string
    :returns: list of fault dictionaries
    :rtype: list
    """
    fault_list = [];
    [lon, lat, depth, strike, dip, length, width, ll_slip,
     thrust_slip, _, _] = np.loadtxt(infile, skiprows=1, unpack=True, dtype={"names": ('lon', 'lat', 'depth', 'strike',
                                                                                       'dip', 'length', 'width', 'ss',
                                                                                       'ds', 'tensile', 'num'),
                                                                             "formats": (np.float, np.float, np.float,
                                                                                         np.float, np.float, np.float,
                                                                                         np.float, np.float, np.float,
                                                                                         np.float, np.float)});
    for i in range(len(lon)):
        one_fault = {"strike": strike[i], "dip": dip[i], "length": length[i] / 1000, "width": width[i] / 1000,
                     "depth": -depth[i] / 1000};
        center_lon = lon[i];
        center_lat = lat[i];
        x_start, y_start = fault_vector_functions.add_vector_to_point(0, 0, one_fault["length"] / 2,
                                                                      one_fault["strike"] - 180);  # in km
        corner_lon, corner_lat = fault_vector_functions.xy2lonlat(x_start, y_start, center_lon, center_lat);
        one_fault["lon"] = corner_lon;
        one_fault["lat"] = corner_lat;
        one_fault["rake"] = fault_vector_functions.get_rake(ll_slip[i], thrust_slip[i]);
        one_fault["slip"] = fault_vector_functions.get_total_slip(ll_slip[i], thrust_slip[i]);
        fault_list.append(one_fault);
    return fault_list;


def write_faults_intxt(faults, outfile, receiver=True, source=False, write_header=False):
    """
    Writes faults to intxt as receivers with zero slip, or sources with finite slip.

    :param faults: list of fault dictionaries
    :type faults: list
    :param outfile: name of output intxt file
    :type outfile: string
    :param receiver: whether the faults are receivers (slip set to 0), default True
    :type receiver: bool
    :param source: whether the faults are sources (finite slip), default False
    :type source: bool
    :param write_header: whether to write the header line, default True
    :type write_header: bool
    """
    ofile = open(outfile, 'w');
    ofile.write(
        '#S/R: strike rake dip length(km) width(km) updip_corner_lon updip_corner_lat updip_corner_dep(km) slip\n')
    if write_header:
        ofile.write("General: 0.250 0.40 -115.75 -115.25 -115.5 32.75 33.25 33.0 # general information\n");
        ofile.write("Receiver: 191 -103 45 6 7 -115.507 33.1 0  # random fault\n");
    for fault in faults:
        if receiver is True:
            ofile.write("Receiver: %.1f nan %.1f %.2f %.2f %.3f %.3f %.2f \n" % (
                fault["strike"], fault["dip"], fault["length"], fault["width"], fault["lon"], fault["lat"],
                fault["depth"]));
        if source is True:
            ofile.write("Source_Patch: %.1f %.1f %.1f %.2f %.2f %.3f %.3f %.2f %.3f \n" % (
                fault["strike"], fault["rake"], fault["dip"], fault["length"], fault["width"], fault["lon"],
                fault["lat"],
                fault["depth"], fault["slip"]));
    ofile.close();
    return;


def write_faults_json(faults_list, outfile):
    """
    Writes faults to json as receivers with zero slip

    :param faults: list of fault dictionaries
    :type faults: list
    :param outfile: name of output json file
    :type outfile: string
    """
    output = {};
    faults = {};
    count = 0;
    for fault in faults_list:
        # Convert the fault (which has top left corner) into a fault with top center coordinate
        corner_lon = fault["lon"]
        corner_lat = fault["lat"]
        x_center, y_center = fault_vector_functions.add_vector_to_point(0, 0, fault["length"] / 2, fault["strike"]);
        center_lon, center_lat = fault_vector_functions.xy2lonlat(x_center, y_center, corner_lon, corner_lat);

        count = count + 1;
        label = "fault" + str(count)
        fault["length"] = fault["length"] * 1000;
        fault["width"] = fault["width"] * 1000;
        fault["basis1"] = [1, 0, 0];
        fault["basis2"] = None;
        fault["Nlength"] = 1;
        fault["Nwidth"] = 1;
        fault["penalty"] = 1;
        fault["position"] = [center_lon, center_lat, fault["depth"]];
        fault.pop("lon");
        fault.pop("lat");
        fault.pop("depth");
        faults[label] = fault;
    output["faults"] = faults;
    output["plotter"] = "gmt";
    with open(outfile, 'w') as ofile:
        json.dump(output, ofile, indent=4);
    return;


def write_slipdistribution(faults, outfile):
    """Write a slip distribution for Slippy.
    lon[deg] lat[deg] depth[m] strike[deg] dip[deg] len[m] wid[m] left-lateral[m] thrust[m] tensile[m] segment_num.
    Depths are negative in meters. lon/lat refer to fault center.

    Function not written yet.
    """
    ofile = open(outfile, 'w');
    ofile.write(
        '# lon[degrees] lat[degrees] depth[m] strike[degrees] dip[degrees] length[m] width[m] left-lateral[m] thrust['
        'm] tensile[m] segment_num\n');
    # for fault in faults:
    print("THIS FUNCTION ISN'T WRITTEN YET");
    ofile.close();
    return;


# ------- CONVERSION FUNCTIONS ------------ # 

def json2intxt(infile, outfile):
    """Converts json faults to .intxt faults"""
    faults = read_faults_json(infile);
    write_faults_intxt(faults, outfile, receiver=True, source=False);
    return;


def intxt2json(infile, outfile):
    """Converts .intxt faults to json faults"""
    faults = read_faults_intxt(infile);
    write_faults_json(faults, outfile);
    return;


def slippydist2intxt(infile, outfile):
    """Converts slippy faults to .intxt faults"""
    faults = read_slippy_distribution(infile);
    write_faults_intxt(faults, outfile, receiver=False, source=True, write_header=True);
    return;


def intxt2slipdistribution(infile, outfile):
    """Converts .intxt faults to slippy faults (in progress)"""
    faults = read_faults_intxt(infile);
    write_slipdistribution(faults, outfile);
    return;

