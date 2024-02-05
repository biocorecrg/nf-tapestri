# nf-tapestri
Nextflow pipeline for Beekiman lab

## Description graph

![GRAPH](https://github.com/biocorecrg/nf-tapestri/blob/main/img/graph_tapestri.png?raw=true)

## Road map
- [x] little test dataset
- [x] small prototype with one step
- [x] adding more steps
- [ ] adding final report
- [ ] writing documentation


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

### Downloading the tapestry image
For downloading the tapestri singularity image you can use the INSTALL.sh script that essentially does a `wget` and a `move` to the path with all the singularity directory.

```
sh INSTALL.sh $user $password $path
```

Un our HPC the cache folder for the singularity images is:

```
/software/bi/biocore_tools/singularity_containers/
```

But you can change it by editing the `nextflow.config` file at the section `singularity.cacheDir`


## Update

For upgrading the pipeline just run inside the `nf-tapestri` folder

```
git pull --recurse-submodules 
```


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
nreads: 100000
hook: ""
email: ""
output: "${projectDir}/outfolder"
```

The command line of the single tools can be tweaked using the file `tools_opt.tsv`. This is the one currently:

```
#tool	extrapars

#tool	extrapars
tapestri-barcode-extract	"--chemistry dna.V2 --no-cut-adapters --decode --barcode-suffix -1"
bwa	"-C -M"
samtools	"-bh -q 30 -F 4 -F 8 -F 0X0100"
picard	"SORT_ORDER=coordinate"
cellfinder_MCL_HC.R	"-c 0.7"

```

## Running the pipeline

```
nextflow run main.nf -params-file params.yaml -bg > log.txt
```
We recommend redirecting the `work` folder to a space without backup like `scratch` (you have to ask SIT for having a folder there).

```
nextflow run main.nf -params-file params.yaml -bg -w //nfs/scratch01/YOURGROUP/work > log.txt
```

## Ouptut folders

The following folders are generated:

- bam
- barcode
- distribution
