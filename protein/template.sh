#!/bin/bash
#$ -cwd 
#$ -q short-centos79,short-rocky9,long-centos79,long-centos79
#$ -pe smp 8
#$ -j y
#$ -l virtual_free=80G
#$ -l h_rt=6:00:00
#source /users/lvelten/mscherer/software/anaconda3/bin/activate /users/lvelten/mscherer/conda/envs/tapestri ;
time singularity exec -e /software/bi/biocore_tools/singularity_containers/biocorecrg-tapestri-2.0.1.img tapestri protein run \
-w 1 \
--ab-adapters-3 ab_adapters_3.fasta \
--ab-adapters-5 ab_adapters_5.fasta \
--ab-barcodes totalseq-d-heme-oncology-CD27.csv \
--output-folder ${3} \
--output-prefix T41 \
--n-cores 8 \
--r1 ${1} \
--r2 ${2}

