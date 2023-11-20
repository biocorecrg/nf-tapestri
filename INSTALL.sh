#!/bin/bash
#
# This script will download the singularity image with tapestri 
# params:
# 1 user
# 2 password

if [ x"$1" == x ]; then
        echo "please enter username"
	exit
fi
if [ x"$2" == x ]; then
        echo "please enter password"
        exit
fi
if [ x"$3" == x ]; then
        echo "please define the path where the store the image"
        exit
fi

wget --user $1 --password $2 https://public-docs.crg.es/biocore/projects/tapestri/biocorecrg-tapestri-2.0.1.img
mv biocorecrg-tapestri-2.0.1.img $3/
