#!/bin/bash
#$ -cwd 
#$ -q short-centos79,short-rocky9,long-centos79,long-centos79
#$ -pe smp 1
#$ -j y
#$ -l virtual_free=20G
#$ -l h_rt=6:00:00
zcat $1 | split -l $3 - $2



