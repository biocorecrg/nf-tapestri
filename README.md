# nf-tapestri
Nextflow pipeline for Beekiman lab

## Description graph

![GRAPH](https://github.com/biocorecrg/nf-tapestri/blob/main/img/graph_tapestri.png?raw=true)

## Road map
1. little test dataset
2. small prototype with one step
3. adding more steps
4. adding final report
5. writing documentation


## Install
- Use singularity
Edit the `$HOME/.bashrc` file adding:

```
module use /software/as/el7.2/EasyBuild/CRG/modules/all
module load Java/11.0.2
module load Singularity/3.7.0
```

- Install nextflow

```
curl -s https://get.nextflow.io | bash
```

- Install the pipeline

```
git clone git@github.com:biocorecrg/nf-tapestri.git --recurse-submodules 
```

or 

```
git clone https://github.com/biocorecrg/nf-tapestri.git --recurse-submodules 
```

## Upgrade

For upgrading the pipeline just run inside the `nf-tapestri` folder

```
git pull --recurse-submodules 
```


## Parameters
Pipeline parameters are inside the `params.yaml` file:

```
# Parameters
reads: "${projectDir}/test_data/*_read{1,2}.fastq.gz"
hook: ""
email: ""
output: "${projectDir}/outfolder"

```

The command line of the single tools can be tweaked using the file `tools_opt.tsv`. This is the one currently:

```
#tool	extrapars
tapestri-barcode-extract	"--chemistry dna.V2 --no-cut-adapters --decode --barcode-suffix -1"

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

...
