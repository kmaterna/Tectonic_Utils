# The functions in this script convert between two fault formats 
# Format 1: json, Trever Hines' format for Slippy
# Format 2: .intxt, Kathryn Materna's format for Elastic_stresses_py
# Fault is a dictionary here containing: strike, dip, length(km), width(km), lon, lat, depth
# Format  json: basis1, basis2, length(m), width(m), nlength, nwidth, strike, dip, position [lon, lat, dep], penalty
# Format intxt: strike rake dip length(km) width(km) updip_corner_lon updip_corner_lat updip_corner_dep
# Lon, lat, depth refer to top left corner of the fualt. 

import numpy as np 
import json
import haversine

def read_faults_intxt(infile):
	# Read all faults that are sources or receivers
	fault_list=[];
	ifile=open(infile,'r');
	for line in ifile:
		temp=line.split();
		if temp[0]=="S:" or temp[0]=="R:":
			one_fault={};
			one_fault["strike"]=float(temp[1]);
			one_fault["dip"]=float(temp[3]);
			one_fault["length"]=float(temp[4]);
			one_fault["width"]=float(temp[5]);
			one_fault["lon"]=float(temp[6]);
			one_fault["lat"]=float(temp[7]);
			one_fault["depth"]=float(temp[8]);
			fault_list.append(one_fault);
	ifile.close();
	return fault_list;


def read_faults_json(infile):
	# Read all faults from a json file (just geometry; don't have slip or rake)
	# It has to convert from fault center to fault corner. 
	fault_list=[];
	config_file = open('config.json','r')
	config = json.load(config_file);
	for key in config["faults"].keys():
		one_fault={};
		one_fault["strike"]=config["faults"][key]["strike"];
		one_fault["dip"]=config["faults"][key]["dip"];
		one_fault["length"]=config["faults"][key]["length"]/1000.0;
		one_fault["width"]=config["faults"][key]["width"]/1000.0;
		center_lon = config["faults"][key]["position"][0];
		center_lat = config["faults"][key]["position"][1];
		x_start, y_start = add_vector_to_point(0, 0, one_fault["length"]/2, -one_fault["strike"]);  # in km
		corner_lon, corner_lat = xy2lonlat(xi,yi,center_lon,center_lat);  # 
		one_fault["lon"]=corner_lon;
		one_fault["lat"]=corner_lat;
		one_fault["depth"]=config["faults"][key]["position"][2];
		fault_list.append(one_fault);
	return fault_list;


def write_faults_intxt(faults, outfile):
	# Writes the files as receivers with zero slip. 
	# Rake is unknown because we only have geometry. 
	ofile=open(outfile,'w');
	for fault in faults:
		ofile.write("R: %.1f nan %.1f %.2f %.2f %.3f %.3f %.2f \n" % (fault["strike"],fault["dip"],fault["length"],fault["width"],fault["lon"],fault["lat"],fault["depth"]) );
	ofile.close();
	return;

def write_faults_json(faults_list, outfile):
	output={};
	faults={};
	count=0;
	for fault in faults_list:
		# Convert the fault (which has top left corner) into a fault with top center coordinate
		corner_lon = fault["lon"]
		corner_lat = fault["lat"]
		x_center, y_center = add_vector_to_point(0, 0, fault["length"]/2, fault["strike"]);  # in km
		center_lon, center_lat = xy2lonlat(x_center,y_center,corner_lon,corner_lat);  # 

		count=count+1;
		label="fault"+str(count)
		fault["length"]=fault["length"]*1000;
		fault["width"]=fault["width"]*1000;
		fault["basis1"]=[1,0,0];
		fault["basis2"]=None;
		fault["Nlength"]=1;
		fault["Nwidth"]=1;
		fault["penalty"]=1;
		fault["position"]=[center_lon,center_lat,fault["depth"]];
		fault.pop("lon");
		fault.pop("lat");
		fault.pop("depth");
		faults[label]=fault;
	output["faults"]=faults;
	output["plotter"]="gmt";
	with open(outfile,'w') as ofile:
		json.dump(output,ofile,indent=4);
	return;

def json2intxt(infile, outfile):
	faults = read_faults_json(infile);
	write_faults_intxt(faults, outfile);
	return;

def intxt2json(infile, outfile):
	faults = read_faults_intxt(infile);
	write_faults_json(faults, outfile);
	return;


def add_vector_to_point(x0,y0,vector_mag,vector_heading):
	# Vector heading defined as strike- CW from north.
	theta=np.deg2rad(90-vector_heading);
	x1 = x0 + vector_mag*np.cos(theta);
	y1 = y0 + vector_mag*np.sin(theta);
	return x1, y1;

def xy2lonlat(xi,yi,reflon,reflat):
	lat=reflat+( yi*1/111.000 );
	lon=reflon+( xi*1/(111.000*abs(np.cos(np.deg2rad(reflat)))) );
	return lon, lat;


def latlon2xy(loni,lati,lon0,lat0):
	# returns the distance between a point and a reference in km. 
	radius = haversine.distance([lat0,lon0], [lati,loni]);
	bearing = haversine.calculate_initial_compass_bearing((lat0, lon0),(lati, loni))
	azimuth = 90 - bearing;
	x = radius * np.cos(np.deg2rad(azimuth));
	y = radius * np.sin(np.deg2rad(azimuth));
	return [x, y];

