outputBarcode   = "${params.output}/Barcode"


params.EXTRAPARS = ""
params.LABEL = ""

// Create a channel for tool options
def getParameters(pars_tools_file) {
	def pars_tools = file(pars_tools_file)
	if( !pars_tools.exists() ) exit 1, "Missing tools options config: '$pars_tools'"

	def progPars = [:]
	def allLines  = pars_tools.readLines()

	for( line : allLines ) {
    	def list = line.split("\t")
    	if (list.length <1) {
			 error "ERROR!!! Tool option file has to be tab separated\n" 
		}
    	if (!(list[0] =~ /#/ )) {
			progPars[list[0]] = list[1].replace("\"", "").replace('$projectDir', "${projectDir}").replace('${projectDir}', "${projectDir}")
    	}  
	}	
	return(progPars)
}

process barcode_extract {
	publishDir outputBarcode, mode: 'copy'
	container 'biocorecrg/tapestri:2.0.1'
    label (params.LABEL)
	tag("${sampleID}")
	
	input:
	tuple val(sampleID), path(reads)

	output:
	tuple val(sampleID), path("${sampleID}_read1.with.barcodes.fastq.gz"), path("${sampleID}_read2.with.barcodes.fastq.gz"), emit: barcode_reads  
	tuple val(sampleID), path("${sampleID}_raw.barcodes.txt"), emit: raw_barcodes  
	tuple val(sampleID), path("${sampleID}_01_barcode.metrics.json"), emit: barcode_metrics  

	script:
	"""
	tapestri barcode extract \
		--input-fwd ${reads[0]} \
		--input-rev ${reads[1]} \
		--output-fwd ${sampleID}_read1.with.barcodes.fastq.gz \
		--output-rev ${sampleID}_read2.with.barcodes.fastq.gz \
		--output-barcodes ${sampleID}_raw.barcodes.txt \
		--output-metrics ${sampleID}_01_barcode.metrics.json \
		--cores ${task.cpus} ${params.EXTRAPARS}
	"""

}


