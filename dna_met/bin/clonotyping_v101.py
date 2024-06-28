#!/usr/bin/env python
# coding: utf-8

# After the data has been processed using Trust4 to identify the v,d,j segments of both heavy and light chains per cell; we want to clonotype the data.
# 
# After multiple tests, we decided to follow this approach:
#     - Set minimum threshold for input in scirpy to 4 reads (i.e. 4 or more), based on the number of reads
#     - Take out all input barcodes reads that have only light chain or only heavy chain information after this filter
#     - Run Scirpy with Hamming distance 2, receptor arms “all”, same_v_gene on and dual_ir “primary only”. 
# 
# Taking out VDJ and VJ orphan barcodes is majorly based on the fact that:
#     - The majority of our cells present as single pairs
#     - When grouping cells based on the same (=Hamming distance 0) light chain CDR3 gives the presence of multiple V fragments at the heavy chain locus
#     - When grouping cells based on the same (=Hamming distance 0) heavy chain CDR3 gives the presence of multiple V fragments at the light chain locus . The latter two mean that based on heavy or light chain alone our data cannot be trusted.

# ##### Import packages & libraries.
# --> Check what need to be removed later !


#!/usr/bin/env python
# coding: utf-8
import scirpy as ir
import scanpy as sc
from glob import glob
import pandas as pd
import tarfile
import anndata
import warnings
import os
import argparse
from matplotlib import pyplot as plt
import collections
import numpy as np
from matplotlib import cm
sc.settings.verbosity = 4 # verbosity level of scanpy


##### Define and set arguments.

parser = argparse.ArgumentParser(description="Script to identify cells' clonotypes based on both light and heavy chains.")
parser.add_argument("-tf","--trust_file",type=str, required=True, help="Output file from trust4 in the 10x CellRanger vdj format format; usually following this pattern <*_b.csv> or <filtered_contig_annotations.csv>.")
parser.add_argument("-rc","--read_cutoff", default=4, type=str, required=True, help="Minimum number of reads required to keep vdj info.") ## to re-check !!!!!
parser.add_argument("-of","--output_folder", type=str, required=True, help="Path to the output folder.")
parser.add_argument("-sn","--sample_name", type=str, required=True, help="Path to the output folder.")
parser.add_argument("-bc","--cell_barcodes", type=str, help="List of cells barcodes to keep per sample.", required=False, default="")

args = parser.parse_args()
trust_file = args.trust_file
read_cutoff = int(args.read_cutoff)
output_folder = args.output_folder
sample = args.sample_name
cell_barcodes = args.cell_barcodes


wdir = "/home/htoukabri/Data/VDJ/" # ur working directory
# list_samples = ["T26_7"] #"P29_2", "T27_1", "T26_7", "T28_3" , "17_1396" "17_1396_vs_P_29_2"

    
################################################################
########## PREPARING & READING INPUT DATA + METADATA ###########
################################################################
print("......................... Processing sample : ", sample, " .........................")

# Load the BCR data.
adata = ir.io.read_10x_vdj(wdir + sample + "/" + trust_file)
#adata.var_names_make_unique()
#adata
# print(adata.obs_keys)
# print(adata.obs_names)

# Read in the metadata (reads information).
metad = pd.read_csv(wdir + sample + "/" + trust_file, header=0, sep=",")

# Adding sample information.
adata.obs["Sample"] = sample

# Create folder to save plots / final tables.
outpath = os.path.join(output_folder, sample, "scirpy_output")
if not os.path.exists(outpath):
    os.makedirs(outpath)

################################################################
#################### GENERAL QC & FILTERING ####################
################################################################
# QC.
ir.tl.chain_qc(adata) #Perform quality control based on the receptor-chain pairing configuration. Categorizes cells into their receptor types and according to their chain pairing status. The function adds three columns to adata.obs, two containing a coarse and fine annotation of receptor types, a third classifying cells according to the number of matched receptor types.
plt.figure(figsize=(12, 9))

ir.pl.group_abundance(adata, groupby="receptor_subtype", target_col="Sample", fig_kws = {'figsize': (6, 4), 'dpi': 120}) # to adjust figsize
plt.savefig(outpath + "/00_receptor_subtype.png", dpi=300)
#plt.show()

print(
    "Fraction of cells with more than one pair of BCRs: {:.2f}".format(np.sum(
        adata.obs["chain_pairing"].isin(["extra VJ", "extra VDJ", "two full chains"]))/ adata.n_obs)
) # Fraction of cells with more than one pair of BCRs: 0.01


# Filtering.
# Remove the ones with multi-chain as they can be doublets.
adata = adata[adata.obs["multi_chain"] != "True", :].copy()

# If a list of barcodes is provided, keep only the given BCs.
if cell_barcodes:
    barcodes = pd.read_csv(wdir + sample + "/" + cell_barcodes)
    print("Number of cells in the sample {} is : {:d}".format(sample, barcodes.shape[0]))
    barcodes.iloc[:,0] = barcodes.iloc[:,0].str.replace("-1", "", regex= True) # remove the trailing "-1" at the end if it exists.
    adata = adata[adata.obs.index.isin(barcodes.iloc[:,0])]

# Keep cells where BCR is detected and filter for doublets : ** maybe not even needed here !**
adata = adata[adata.obs['receptor_type'] == "BCR"]

# Keep only the ones with min reads n == cutoff.
metad = metad[metad["reads"]>= read_cutoff] 
# calculate the sum of reads per cell supporting the vdj info.
reads_pbc = metad.groupby(by="barcode")["reads"].sum() # barcode is already set as .index
#reads_pbc.index = reads_pbc.barcode

# len(reads_pbc[reads_pbc>= 8 ]) #9856
# len(reads_pbc[reads_pbc>= 4 ]) #11483
# >> we'll keep the barcodes with min 8 reads; with 4 for each chain !!!! 

reads_pbc = reads_pbc[reads_pbc>= (read_cutoff*2) ]
# barcodes4reads = metad[metad["reads"]>= read_cutoff]["barcode"].unique()    
barcodes4reads = reads_pbc.index.unique()
adata = adata[adata.obs.index.isin(barcodes4reads)]    

## Optiannaly we can add the reads info to the final table (or just use the BCs that pass the threshhold).
adata.obs["reads_p_bc"] = reads_pbc

# Keep only cells with a pair of chains (Heavy & light).
print("#Number per cells within category: ", adata.obs.groupby(by="chain_pairing")["chain_pairing"].count()) # Just to chack the number of cells with each category.
adata= adata[adata.obs["chain_pairing"] == "single pair"]

################################################################
################## CLONOTYPE DEFINITION ########################
################################################################
## For more details about scirpy parameters:
### https://scirpy.scverse.org/en/latest/generated/scirpy.tl.define_clonotypes.html


# Set parameters for the clonotype definition. !!(We can make it part of arguments, if needed)!!
dist = "hamming"
cutoff = 2
receptor = "all"
dual_chains = "primary_only"
v_gene = True
    
# Compute a sequence-distance metric between all unique VJ CDR3 sequences and between all unique VDJ CDR3 sequences. Default metric = "identity" & Default sequence= "nt"
ir.pp.ir_dist(adata, sequence="nt", metric= dist, cutoff= cutoff, n_jobs=4)  # will be stored in "ir_dist_nt_hamming"

# Cluster clonotypes based on the computed distance.
ir.tl.define_clonotype_clusters(adata, sequence="nt", metric=dist, 
            receptor_arms= receptor, dual_ir=dual_chains, same_v_gene=v_gene, within_group='receptor_type', 
            partitions='connected', resolution=1, n_iterations=5, inplace=True, n_jobs=4) # created key : "cc_nt_hamming"

# Plot the number of cells per clonotype (first 15 clones).
# adata.obs.rename(columns={"cc_nt_hamming": "clone_id", "cc_nt_hamming_size": "clone_size"}, inplace=True)
# sorted_sizes = adata.obs.groupby(by="cc_nt_hamming")["cc_nt_hamming_size"].sort_values(ascending=False)
sorted_sizes = adata.obs[["cc_nt_hamming", "cc_nt_hamming_size"]].drop_duplicates().sort_values(by= "cc_nt_hamming_size", ascending=False)
top_15_sizes = sorted_sizes.iloc[:15]
# top_15_sizes.index = top_15_sizes["cc_nt_hamming"]

fig = plt.figure(figsize = (10, 7))
plt.bar(top_15_sizes["cc_nt_hamming"], top_15_sizes["cc_nt_hamming_size"], color="teal", width=0.4)
plt.xlabel("Clone ID")
plt.xticks(rotation=30)
plt.ylabel("Clone size")
plt.title("Top 15 clones")
plt.savefig(outpath + "/01_clone_sizes.png", dpi=300)
# plt.savefig("/home/htoukabri/Data/VDJ/01_474642_clone_sizes.png", dpi=300)

################################################################
#################### SAVING OUTPUT DATA ########################
################################################################

adata.obs.to_csv(outpath + "/scirpy_obs_out.csv", sep = "\t")

## To run the command; u can do : python3 clonotyping_v101.py -tf "_b.csv" -rc 4 -of . -sn "T26_7"
