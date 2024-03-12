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


if [ x"$3" = x ]; then

        echo "please specify the number of threads"

        exit 1

fi

samtools view -@ $3 $1 | grep -oP 'RG:Z:\K.+?\b' | sort | uniq -c | sort -n -r | awk '{{$1=$1}}1' > $2


 
