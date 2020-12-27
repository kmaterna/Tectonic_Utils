# Netcdf reading and writing functions
# For the InSAR library, only Netcdf3 and Netcdf4 files with PIXEL NODE REGISTRATION are valid.
# The assumption is 2D Netcdf files with 3 variables, in x-y-z order.

import numpy as np
import scipy.io.netcdf as netcdf
import datetime as dt
import subprocess, sys
from netCDF4 import Dataset


# --------------- READING ------------------- #

def parse_pixelnode_registration(filename):
    output = subprocess.check_output(['gmt', 'grdinfo', filename], shell=False);
    if "Pixel node registration used" not in str(output):
        print("ERROR! %s is not a pixel-node registered grd file (pixel-node required for this library)" % filename);
        sys.exit(1);
    return;


def read_netcdf3(filename):
    # A general netcdf3 function that may or may not take variables.
    # Spits out the right thing based on several known patterns (x, y, z; lon,lat,z; etc)
    # Imposes pixel-node registration.
    print("Reading file %s " % filename);
    file = netcdf.netcdf_file(filename, 'r');
    parse_pixelnode_registration(filename);
    [xkey, ykey, zkey] = file.variables.keys();
    xdata0 = file.variables[xkey][:];
    ydata0 = file.variables[ykey][:];
    zdata0 = file.variables[zkey][::];
    return [xdata0.copy(), ydata0.copy(), zdata0.copy()];


def read_netcdf4(filename):
    # Reading a generalized netCDF4 file with 3 variables, using netCDF4 library
    # Examples: (x,y,z; lon,lat,z; etc.)
    # Should read such that the data point is the center of each pixel
    print("Reading file %s " % filename);
    rootgrp = Dataset(filename, "r");
    parse_pixelnode_registration(filename);
    [xkey, ykey, zkey] = rootgrp.variables.keys()
    xvar = rootgrp.variables[xkey];
    yvar = rootgrp.variables[ykey];
    zvar = rootgrp.variables[zkey];
    return [xvar[:], yvar[:], zvar[:, :]];


def read_any_grd(filename):
    # Switch between netcdf4 and netcdf3 automatically.
    try:
        [xdata, ydata, zdata] = read_netcdf3(filename);
    except TypeError:
        [xdata, ydata, zdata] = read_netcdf4(filename);
    return [xdata, ydata, zdata];


def give_metrics_on_grd(filename):
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
    filereader = netcdf.netcdf_file(filename, 'r');
    tdata0 = filereader.variables['t'][:];
    xdata0 = filereader.variables['x'][:]
    ydata0 = filereader.variables['y'][:];
    zdata0 = filereader.variables['z'][:, :, :];
    return [tdata0.copy(), xdata0.copy(), ydata0.copy(), zdata0.copy()];


# --------------- WRITING ------------------- # 


def write_temp_output_txt(z, outfile):
    # A helper function for dumping grid data into pixel-node-registered grd files
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
    # Writing PIXEL NODE registered netcdf4 file from numpy array
    # strategy: send out to a file and make GMT convert to netcdf
    print("writing outfile %s " % outfile);
    outtxt = outfile+'.xyz'
    write_temp_output_txt(z, outtxt);
    xinc = x[1] - x[0];
    yinc = y[1] - y[0];
    xmin = np.min(x)-xinc/2;  # for the half-pixel outside the edge
    xmax = np.max(x)+xinc/2;  # required when writing pixel-node registration from Python's netcdf into .grd files
    ymin = np.min(y)-yinc/2;  # required when writing pixel-node registration from Python's netcdf into .grd files
    ymax = np.max(y)+yinc/2;  # required when writing pixel-node registration from Python's netcdf into .grd files
    increments = str(xinc)+'/'+str(yinc);
    region = str(xmin)+'/'+str(xmax)+'/'+str(ymin)+'/'+str(ymax);
    command = 'gmt xyz2grd '+outtxt+' -G'+outfile+' -I'+increments+' -R'+region+' -ZBLf -r -fg -di-9999 '
    print(command);
    subprocess.call(['gmt', 'xyz2grd', outtxt, '-G'+outfile, '-I'+increments, '-R'+region, '-ZBLf', '-r',
                     '-fg', '-di-9999'], shell=False);
    subprocess.call(['rm', outtxt], shell=False);
    return;


def produce_output_netcdf(xdata, ydata, zdata, zunits, netcdfname, dtype=float):
    # Write netcdf3 grid file.
    # NOTE: The pixel vs gridline registration of this function is not guaranteed.
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
    # IF WE NEED TO FLIP DATA:
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
    # Ultimately we will need a function that writes a large 3D array.
    # Each 2D slice is the displacement at a particular time, associated with a time series.
    # zdata comes in as a 2D array where each element is a timeseries (1D array).
    # It must be re-packaged into a 3D array before we save it.
    # Broke during long SoCal experiment for some reason. f.close() didn't work.

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
