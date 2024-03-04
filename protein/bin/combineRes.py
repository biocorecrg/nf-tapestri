#!/usr/bin/env python 
__author__ = 'luca.cozzuto@crg.eu'
# -*- coding utf-8 -*-

#MODULES
import sys
import re
import optparse
from collections import defaultdict
import pprint
import copy 
import os
import glob 


#BODY FUNTIONS
def options_arg():
	usage = "usage: %prog -f <file list > -o <output table>"
	parser = optparse.OptionParser(usage=usage)
	parser.add_option('-f', '--files', help='Input file list', dest="input" )
	parser.add_option('-o', '--output',help='ouput table file', dest="wotus" )

	(opts,args) = parser.parse_args()
	if opts.input and opts.wotus:pass
	else: parser.print_help()
	return (opts)
def __main__ ():
	parseFiles(opts.input, opts.wotus)

#AUXILIAR MODULES
def parseFiles(filelist, ofile, ):
	pp = pprint.PrettyPrinter(indent=4)
	outres = defaultdict(dict)

	header = ""
	if(re.search(',',filelist)): 
		ifiles = filelist.split(",")
	else:
		ifiles = glob.glob(filelist, recursive=False)
	# prepare files files
	for ifile in ifiles:
		nrow = 0
		infile = open(ifile, 'r')
		for line in infile:
			line = line.rstrip()
			nrow = nrow + 1
			if (nrow == 1):
				header = line
			else:
				fields = line.split('\t')
				cell_ab = fields[0] + "---" + fields[1]	
				raw_c = int(fields[2])
				if (cell_ab not in outres):
					outres[cell_ab] = raw_c
				else:
					outres[cell_ab] = outres[cell_ab] + raw_c
		infile.close()
	
	outstr = header + "\n" 
	
	for cell_abs, value in outres.items():
		cell_vals = cell_abs.split("---")
		outstr += cell_vals[0] + "\t" + cell_vals[1] + "\t" + str(value) + "\n"
			
	f = open(ofile,'w')
	f.write(outstr)
	f.close()
	
	return


#Calling
opts = options_arg()
if opts.input and opts.wotus:
	__main__()

