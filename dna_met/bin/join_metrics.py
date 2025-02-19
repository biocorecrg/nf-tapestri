#!/usr/bin/env python

import os
import re
import argparse
import json


################################
# Setting an argument parser . #
################################

parser = argparse.ArgumentParser(description="Script join barcode metrics")
parser.add_argument('-jsl','--listOfJson', nargs='+', help="List of json to summarize", required=True)
parser.add_argument("-of","--outputName", type=str, help="Prefix name for the output file")
args = parser.parse_args()

# Retrieve the argument values
list_files = args.listOfJson
outname = args.outputName
outfile = outname + "_bc_metrics.txt"
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


for json_file in list_files:
	jfile = open(json_file,"r")
	data = json.load(jfile)
	for name_res, this_val_res in data.items():
		curr_val_res = my_res[name_res]
		my_res[name_res] = this_val_res + curr_val_res

ofile = open(outfile, "w") 
for name_res, val_res in my_res.items():
	ofile.write(outname + "\t" + name_res + "\t" + str(val_res) + "\n")

print("----------------------------------- && done && -----------------------------------")
