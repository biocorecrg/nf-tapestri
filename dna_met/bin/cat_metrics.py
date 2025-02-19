#!/usr/bin/env python

import os
import re
import argparse
import json


################################
# Setting an argument parser . #
################################

parser = argparse.ArgumentParser(description="Script join barcode metrics")
parser.add_argument('-tsl','--listOftxt', nargs='+', help="List of txt to summarize", required=True)
parser.add_argument("-of","--outputName", type=str, help="Name for the output file")
args = parser.parse_args()

# Retrieve the argument values
list_files = args.listOftxt
outfile = args.outputName
################################
############ FUNCs . ###########
################################

my_res = {
	"n_read_pairs": 0,
	"n_read_pairs_trimmed": 0,
	"n_read_pairs_valid_cell_barcodes": 0,
	"n_bases_r1": 0,
	"n_bases_r2": 0,
	"n_bases_r1_q30": 0,
	"n_bases_r2_q30": 0,
	"n_cell_barcode_bases": 0,
	"n_cell_barcode_bases_q30": 0
}

wfile = open(outfile, "w")
wfile.write("# plot_type: 'table'\n")
wfile.write("# section_name: 'Read Statistics'\n")
wfile.write("# description: 'QC Metrics for Read Processing (in Millions)'\n")
wfile.write("sample " + " ".join(my_res.keys()) + "\n")

for txt_file in list_files:
	with open(txt_file, "r") as file:
		data_dict = {}
		lines = file.readlines()
		for line in lines:
			parts = line.strip().split("\t")  # Tab-separated
			sample, metric, value = parts
			mval = round(int(value)/1000000)
			data_dict[metric] = str(mval)
		wfile.write(sample + " ")
		for res in my_res.keys():
			wfile.write(data_dict[res] + " ")
		wfile.write("\n")

wfile.close()

print("----------------------------------- && done && -----------------------------------")
