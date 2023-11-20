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
	tuple val(sampleID), path("${sampleID}.sam.header"), emit: sam_header

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
	
		get_header.sh ${sampleID}_raw.barcodes.txt ${sampleID}.sam.header
	"""

}

process make_minus_5panel {
	//publishDir outputBarcode, mode: 'copy'
	container 'biocorecrg/tapestri:2.0.1'
    label (params.LABEL)
	tag("${panelfile}")

	input:
	path(panelfile)
	
	output:
	path("${panelfile}_minus_5.bed")
	
	script:
	"""
	awk '{print \$1 "\t" (\$2 + 5) "\t" (\$3 - 5) "\t" \$4}' ${panelfile} > ${panelfile}_minus_5.bed 
	"""
}	

process make_target_counts {
	container 'biocorecrg/tapestri:2.0.1'
    label (params.LABEL)
	tag("${sampleID}")
	
	input:
	tuple val(sampleID), path(reads), path(index)

	output:
	tuple val(sampleID), path("${sampleID}.mapped.target.count.txt")

	script:
	"""
	make_target_counts.sh ${reads} ${sampleID}.mapped.target.count.txt ${task.cpus}
 	"""

}


process make_distribution {
	publishDir outputBarcode, mode: 'copy'
	container 'biocorecrg/tapestri:2.0.1'
        label (params.LABEL)
 	tag("${sampleID}")
	
	input:
	tuple val(sampleID), path(reads), path(index), path(panel5)

	output:
	tuple val(sampleID), path("${sampleID}.cellfinder.barcode.distribution.txt"), emit: dist_read1 
        tuple val(sampleID), path("${sampleID}.cellfinder.barcode.distribution_r2.txt"), emit: dist_read2

	script:
	"""
		/home/biodocker/tapestripipeline/bin/resources/Perl/distribution_barcodes_amplicon.pl \
		--amplicon-input ${panel5} \
		--bam-input ${reads} \
		--output-r1 ${sampleID}.cellfinder.barcode.distribution.txt \
		--output-r2 ${sampleID}.cellfinder.barcode.distribution_r2.txt 
 	"""

}



process cellfinder_MCL {
        publishDir outputBarcode, mode: 'copy'
        container 'biocorecrg/tapestri_r:0.0.1'
        label (params.LABEL)
        tag("${sampleID}")

        input:
        tuple val(sampleID), path(input), path(amplifile)

        output:
        path("output")

        script:
        """
        cellfinder_MCL_HC.R -f ${input} -a ${amplifile} ${params.EXTRAPARS} -o ./output
        """

}



