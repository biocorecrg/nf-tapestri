# nf-tapestri
Nextflow pipeline to generate scTAM-seq output files.

nf-tapestri involves three pipelines: DNAmet, protein and somatic mutations calling pipelines. 
DNAmet and protein pipelines can run independently, while somatic mutations pipeline depends on outfiles ganerates by the other two pipelines.

## Install

- Use singularity
Edit the `$HOME/.bashrc` file adding:

```
module use /software/as/el7.2/EasyBuild/CRG/modules/all
module load Java/11.0.2
module load Singularity/3.7.0
```

login and logout. From now on you can use singularity. 

- Install nextflow

```
curl -s https://get.nextflow.io | bash
```

You can place the executable in a place with other binaries and add that folder to your `$PATH`

- Install the pipeline

```
git clone git@github.com:biocorecrg/nf-tapestri.git --recurse-submodules 
```

or 

```
git clone https://github.com/biocorecrg/nf-tapestri.git --recurse-submodules 
```

## Downloading the tapestry image
For downloading the tapestri singularity image you can use the INSTALL.sh script that essentially does a `wget` and a `move` to the path with all the singularity directory.

```
sh INSTALL.sh $user $password $path
```

Un our HPC the cache folder for the singularity images is:

```
/users/bi/public-docs/projects/tapestri/biocorecrg-tapestri-2.0.1.img
```

But you can change it by editing the `nextflow.config` file at the section `singularity.cacheDir`


## Update

For upgrading the pipeline just run inside the `nf-tapestri` folder

```
git pull --recurse-submodules 
```

# DNA methylation pipeline
This pipeline generates two output files per sample: 
- DNA methylation reads counts matrix
- Cell metadata file with doublets information

## Parameters
Pipeline parameters are inside the `params.yaml` file:

```
# Parameters
reads: "${projectDir}/test_data/*_read{1,2}.fastq.gz"
genome: "${projectDir}/anno/GRCh37.p13.genome.fa.gz"
panel: "${projectDir}/anno/panel_info.bed" 
# ampli can be "" in case of cellfinder_scTAMseq.R and it will take the 4th column of the panel file
#ampli: "${projectDir}/anno/Control_amplicons_cell_lines.csv"
ampli: ""
# cellfinder can be cellfinder_scTAMseq.R / cellfinder_MCL_HC.R 
cellfinder: "cellfinder_scTAMseq.R" 
nreads: 50000000
hook: ""
email: ""
output: "${projectDir}/outfolder"
```
nreads parameter: 
- 1000M files or bigger we recommend setting nreads to 50M reads.
- Smaller files, ~100M reads, 20M reads could be more convenient.   


The command line of the single tools can be tweaked using the file `tools_opt.tsv`. This is the one currently:

```
#tool	extrapars

#tool	extrapars
tapestri-barcode-extract	"--chemistry dna.V2 --no-cut-adapters --decode --barcode-suffix -1"
bwa	"-C -M"
samtools	"-bh -q 30 -F 4 -F 8 -F 0X0100"
picard	"SORT_ORDER=coordinate"
cellfinder_MCL_HC.R	"-c 0.7"
cellfinder_scTAMseq.R "-c 0.5"

```

## Running the pipeline
IMPORTANT: follow this instruction to the run the pipeline at the CRG cluster

1) Login to general nextflow node 

```
ssh username@nextflow.linux.crg.es
```
2) Access to nextflow3 through shh

```
ssh nextflow3.linux.crg.es
```

Command to run the pipeline:

```
nextflow run main.nf -params-file params.yaml -bg > log.txt
```
We recommend redirecting the `work` folder to a space without backup like `scratch` (you have to ask SIT for having a folder there).

```
nextflow run main.nf -params-file params.yaml -bg -w //nfs/scratch01/YOURGROUP/work > log.txt
```

## Output folders

The following folders are generated per sample:

- output
- intermediate
- qc_log






# Protein pipeline
This pipeline generates an output file listing number of reads per antibody per cell

## Parameters
Pipeline parameters are inside the `params.yaml` file:

```
# Parameters
reads: "${projectDir}/small_data/*_R{1,2}_001.fastq.gz"
adapt3: "${projectDir}/anno/ab_adapters_3.fasta"
adapt5: "${projectDir}/anno/ab_adapters_5.fasta"
ab_barcode: "${projectDir}/anno/totalseq-d-heme-oncology-CD27.csv"
nreads: 50000000
hook: ""
email: ""
output: "${projectDir}/outfolder"
```
 
The command line of the single tools can be tweaked using the file `tools_opt.tsv`.


