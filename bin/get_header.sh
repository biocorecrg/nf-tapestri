#!/bin/bash
#
# This script makes the sam header from barcodes extracted by tapestri pipeline

if [ x"$1" == x ]; then

        echo "please specify the name of the input file"

        exit 1

fi


if [ x"$2" = x ]; then

        echo "please specify the name of the output file"

        exit 1

fi

sed -n 's/^\([^\t]*\).*/@RG\tID:\1\tSM:\1/p' $1 > $2 
 
