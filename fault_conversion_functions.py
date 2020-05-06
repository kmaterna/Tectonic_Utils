# The functions in this script convert between two fault formats 
# Format 1: json, Trever Hines' format for Slippy
# Format 2: .intxt, Kathryn Materna's format for Elastic_stresses_py
# Fault is a dictionary here containing: strike, dip, length(km), width(km), lon, lat, depth
# Format  json: basis1, basis2, length(m), width(m), nlength, nwidth, strike, dip, position [lon, lat, dep], penalty
# Format intxt: strike rake dip length(km) width(km) updip_corner_lon updip_corner_lat updip_corner_dep

import numpy as np 
import json

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
	fault_list=[];
	config_file = open('config.json','r')
	config = json.load(config_file);
	for key in config["faults"].keys():
		one_fault={};
		one_fault["strike"]=config["faults"][key]["strike"];
		one_fault["dip"]=config["faults"][key]["dip"];
		one_fault["length"]=config["faults"][key]["length"]/1000.0;
		one_fault["width"]=config["faults"][key]["width"]/1000.0;
		one_fault["lon"]=config["faults"][key]["position"][0];
		one_fault["lat"]=config["faults"][key]["position"][1];
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
		count=count+1;
		label="fault"+str(count)
		fault["length"]=fault["length"]*1000;
		fault["width"]=fault["width"]*1000;
		fault["basis1"]=[1,0,0];
		fault["basis2"]=None;
		fault["nlength"]=1;
		fault["nwidth"]=1;
		fault["penalty"]=1;
		fault["position"]=[fault["lon"],fault["lat"],fault["depth"]];
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
