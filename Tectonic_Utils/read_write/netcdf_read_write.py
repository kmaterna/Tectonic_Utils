"""
Netcdf reading and writing functions.
Only Netcdf3 and Netcdf4 files with PIXEL NODE REGISTRATION are valid.
The assumption is 2D Netcdf files with 3 variables, in x-y-z order.
"""

import numpy as np
import scipy.io.netcdf as netcdf
import datetime as dt
import subprocess
from netCDF4 import Dataset


# --------------- READING ------------------- #

def parse_pixelnode_registration(filename):
    """Ensure pixel node registration for netcdf file"""
    output = subprocess.check_output(['gmt', 'grdinfo', filename], shell=False)
    assert("Pixel node registration used" in str(output)), ValueError("ERROR! "+filename+" not pixel-node registered");
    return;


def properly_parse_variables(key1, key2, key3):
    """
    Set proper ordering for known keys in a netcdf file. Options: [x, y, z]; [lon, lat, z]; [longitude, latitude, z].

    :param key1: names of netcdf variable key
    :type key1: string
    :param key2: names of netcdf variable key
    :type key2: string
    :param key3: names of netcdf variable key
    :type key3: string
    :returns: ordered keys
    :rtype: list
    """
    key_list = [key1, key2, key3];
    # get the x key:
    if 'x' in key_list:
        xkey = 'x';
    elif 'lon' in key_list:
        xkey = 'lon';
    elif 'longitude' in key_list:
        xkey = 'longitude';
    else:
        raise Exception("Xkey not determined");
    # get the y key
    if 'y' in key_list:
        ykey = 'y';
    elif 'lat' in key_list:
        ykey = 'lat';
    elif 'latitude' in key_list:
        ykey = 'latitude';
    else:
        raise Exception("Ykey not determined");
    key_list.pop(key_list.index(xkey));
    key_list.pop(key_list.index(ykey));
    if len(key_list) == 1:
        zkey = key_list[0];
    else:
        raise Exception("Zkey not determined");
    return [xkey, ykey, zkey];


def read_netcdf3(filename):
    """
    A netcdf3 reading function for pixel-node registered files with recognized key patterns.

    :param filename: name of netcdf3 file
    :type filename: string
    :returns: [xdata, ydata, zdata]
    :rtype: list of 3 np.ndarrays
    """
    print("Reading file %s " % filename);
    file = netcdf.netcdf_file(filename, 'r');
    parse_pixelnode_registration(filename);
    [xkey, ykey, zkey] = file.variables.keys();
    [xkey, ykey, zkey] = properly_parse_variables(xkey, ykey, zkey);
    xdata0 = file.variables[xkey][:];
    ydata0 = file.variables[ykey][:];
    zdata0 = file.variables[zkey][::];
    return [xdata0.copy(), ydata0.copy(), zdata0.copy()];


def read_netcdf4(filename):
    """
    A netcdf4 reading function for pixel-node registered files with recognized key patterns.

    :param filename: name of netcdf4 file
    :type filename: string
    :returns: [xdata, ydata, zdata]
    :rtype: list of 3 np.ndarrays
    """
    print("Reading file %s " % filename);
    rootgrp = Dataset(filename, "r");
    parse_pixelnode_registration(filename);
    [xkey, ykey, zkey] = rootgrp.variables.keys()
    [xkey, ykey, zkey] = properly_parse_variables(xkey, ykey, zkey);
    xvar = rootgrp.variables[xkey];
    yvar = rootgrp.variables[ykey];
    zvar = rootgrp.variables[zkey];
    return [xvar[:], yvar[:], zvar[:, :]];


def read_any_grd(filename):
    """
    A general netcdf4/netcdf3 reading function for pixel-node registered files with recognized key patterns.

    :param filename: name of file
    :type filename: string
    :returns: [xdata, ydata, zdata]
    :rtype: list of 3 np.ndarrays
    """
    try:
        [xdata, ydata, zdata] = read_netcdf3(filename);
    except TypeError:
        [xdata, ydata, zdata] = read_netcdf4(filename);
    return [xdata, ydata, zdata];


def give_metrics_on_grd(filename):
    """Print shape, min/max, and NaN metrics on a netcdf grid file"""
    grid_data = read_any_grd(filename)[2];
    nan_pixels = np.count_nonzero(np.isnan(grid_data));
    total_pixels = np.shape(grid_data)[0] * np.shape(grid_data)[1];
    print("Shape of %s is [%d, %d]" % (filename, np.shape(grid_data)[0], np.shape(grid_data)[1]));
    print("Min data is %f " % (np.nanmin(grid_data)));
    print("Max data is %f " % (np.nanmax(grid_data)));
    print(
        "Nans: %d of %d pixels are nans (%.3f percent)" % (nan_pixels, total_pixels, nan_pixels / total_pixels * 100));
    return;


def read_3D_netcdf(filename):
    """
    Reading function for 3D netcdf pixel-node registered files with key pattern 't, x, y, z'

    :param filename: name of netcdf file
    :type filename: string
    :returns: [tdata, xdata, ydata, zdata]
    :rtype: list of 4 np.ndarrays
    """
    filereader = netcdf.netcdf_file(filename, 'r');
    tdata0 = filereader.variables['t'][:];
    xdata0 = filereader.variables['x'][:]
    ydata0 = filereader.variables['y'][:];
    zdata0 = filereader.variables['z'][:, :, :];
    return [tdata0.copy(), xdata0.copy(), ydata0.copy(), zdata0.copy()];


# --------------- WRITING ------------------- # 


def write_temp_output_txt(z, outfile):
    """A helper function for dumping grid data into pixel-node-registered grd files"""
    (y, x) = np.shape(z);
    z = np.reshape(z, (x*y,));
    # z = np.array(z).astype(str)
    # z[z == 'nan'] = '-9999'
    # np.savetxt(outfile, z, fmt='%s');
    ofile = open(outfile, 'w+b');   # binary file
    zbytes = bytearray(z)
    ofile.write(zbytes)
    ofile.close()
    print("writing temporary outfile %s " % outfile);
    return;


def write_netcdf4(x, y, z, outfile):
    """
    Writing PIXEL NODE registered netcdf4 file from numpy array.
    Internal strategy: send out to a binary file and make GMT convert to netcdf.
    """
    print("writing outfile %s " % outfile);
    outtxt = outfile+'.xyz'
    write_temp_output_txt(z, outtxt);
    xinc = np.round(x[1] - x[0], 6);
    yinc = np.round(y[1] - y[0], 6);
    xmin = np.round(np.min(x)-xinc/2, 6);  # for the half-pixel outside the edge
    xmax = np.round(np.max(x)+xinc/2, 6);  # writing pixel-node registration from Python's netcdf into .grd files
    ymin = np.round(np.min(y)-yinc/2, 6);  # writing pixel-node registration from Python's netcdf into .grd files
    ymax = np.round(np.max(y)+yinc/2, 6);  # writing pixel-node registration from Python's netcdf into .grd files
    increments = str(xinc)+'/'+str(yinc);
    region = str(xmin)+'/'+str(xmax)+'/'+str(ymin)+'/'+str(ymax);
    if isinstance(z[0][0], np.float64):
        binary_format_flags = '-ZBLd';   # double precision floating piont number, standard numpy float
    else:
        binary_format_flags = '-ZBLf';   # 4-byte floating point number
    command = 'gmt xyz2grd '+outtxt+' -G'+outfile+' -I'+increments+' -R'+region+' '+binary_format_flags+\
              ' -r -fg -di-9999 '
    print(command);
    subprocess.call(['gmt', 'xyz2grd', outtxt, '-G'+outfile, '-I'+increments, '-R'+region, binary_format_flags, '-r',
                     '-fg', '-di-9999'], shell=False);
    subprocess.call(['rm', outtxt], shell=False);
    return;


def produce_output_netcdf(xdata, ydata, zdata, zunits, netcdfname, dtype=float):
    """
    Write netcdf3 grid file.
    NOTE: The pixel vs gridline registration of this function is not guaranteed;
    depends on file system and float type and precision :(.
    Safer to use write_netcdf4().
    """
    print("Writing output netcdf to file %s " % netcdfname);
    f = netcdf.netcdf_file(netcdfname, 'w');
    f.history = 'Created for a test';
    f.createDimension('x', len(xdata));
    f.createDimension('y', len(ydata));
    print(np.shape(zdata));
    x = f.createVariable('x', dtype, ('x',))
    x[:] = xdata;
    x.units = 'range';
    y = f.createVariable('y', dtype, ('y',))
    y[:] = ydata;
    y.units = 'azimuth';
    z = f.createVariable('z', dtype, ('y', 'x',));
    z[:, :] = zdata;
    z.units = zunits;
    f.close();
    flip_if_necessary(netcdfname);
    return;


def flip_if_necessary(filename):
    """If netcdf3 file is stored with xinc or yinc backwards, we replace with a copy that flips the affected axis."""
    xinc = subprocess.check_output('gmt grdinfo -M -C ' + filename + ' | awk \'{print $8}\'',
                                   shell=True);  # x-increment
    yinc = subprocess.check_output('gmt grdinfo -M -C ' + filename + ' | awk \'{print $9}\'',
                                   shell=True);  # y-increment
    xinc = float(xinc.split()[0]);
    yinc = float(yinc.split()[0]);

    if xinc < 0:  # FLIP THE X-AXIS
        print("flipping the x-axis");
        [xdata, ydata] = read_netcdf3(filename)[0:2];
        data = read_netcdf3(filename)[2];
        # This is the key! Flip the x-axis when necessary.
        # xdata=np.flip(xdata,0);  # This is sometimes necessary and sometimes not!  Not sure why.
        produce_output_netcdf(xdata, ydata, data, 'mm/yr', filename);
        xinc = subprocess.check_output('gmt grdinfo -M -C ' + filename + ' | awk \'{print $8}\'',
                                       shell=True);  # x-increment
        xinc = float(xinc.split()[0]);
        print("New xinc is: %f " % xinc);
    if yinc < 0:
        print("flipping the y-axis");
        [xdata, ydata] = read_netcdf3(filename)[0:2];
        data = read_netcdf3(filename)[2];
        # Flip the y-axis when necessary.
        # ydata=np.flip(ydata,0);
        produce_output_netcdf(xdata, ydata, data, 'mm/yr', filename);
        yinc = subprocess.check_output('gmt grdinfo -M -C ' + filename + ' | awk \'{print $9}\'',
                                       shell=True);  # y-increment
        yinc = float(yinc.split()[0]);
        print("New yinc is: %f" % yinc);
    return;


def produce_output_TS_grids(xdata, ydata, zdata, timearray, zunits, outdir):
    """ Write many netcdf3 files, one for each step of a timearray. Each file will be named with a datetime suffix."""
    print("Shape of zdata originally:", np.shape(zdata));
    for i in range(len(timearray)):
        filename = dt.datetime.strftime(timearray[i], "%Y%m%d") + ".grd";
        zdata_slice = np.zeros([len(ydata), len(xdata)]);
        for k in range(len(xdata)):
            for j in range(len(ydata)):
                temp_array = zdata[j][k][0];
                zdata_slice[j][k] = temp_array[i];
        produce_output_netcdf(xdata, ydata, zdata_slice, zunits, outdir + "/" + filename);
    return;


def produce_output_timeseries(xdata, ydata, zdata, timearray, zunits, netcdfname):
    """"
    Write dataset with t, x, y, z into large 3D netcdf.
    Each 2D slice is the displacement at a particular time, associated with a time series.
    zdata comes in as a 2D array where each element is a timeseries (1D array), so it must be re-packaged into
    3D array before we save it.
    Broke during long SoCal experiment for some reason. f.close() didn't work.
    """

    print("Shape of zdata originally:", np.shape(zdata));
    zdata_repacked = np.zeros([len(timearray), len(ydata), len(xdata)]);
    print("Intended repackaged zdata of shape: ", np.shape(zdata_repacked));
    if np.shape(zdata) == np.shape(zdata_repacked):
        print("No repacking necessary");
        zdata_repacked = zdata;
    else:
        print("Repacking zdata into zdata_repacked");
        for i in range(len(zdata[0][0][0])):  # for each time interval:
            print(i);
            for k in range(len(xdata)):
                for j in range(len(ydata)):
                    temp_array = zdata[j][k][0];
                    zdata_repacked[i][j][k] = temp_array[i];

    print("Writing output netcdf to file %s " % netcdfname);
    days_array = [];
    for i in range(len(timearray)):
        delta = timearray[i] - timearray[0];
        days_array.append(delta.days);
    f = netcdf.netcdf_file(netcdfname, 'w');
    f.history = 'Created for a test';
    f.createDimension('t', len(timearray));
    f.createDimension('x', len(xdata));
    f.createDimension('y', len(ydata));

    t = f.createVariable('t', 'i4', ('t',))
    t[:] = days_array;
    t.units = 'days';
    x = f.createVariable('x', float, ('x',))
    x[:] = xdata;
    x.units = 'range';
    y = f.createVariable('y', float, ('y',))
    y[:] = ydata;
    y.units = 'azimuth';

    z = f.createVariable('z', float, ('t', 'y', 'x'));
    z[:, :, :] = zdata_repacked;
    z.units = zunits;
    f.close();
    return;
