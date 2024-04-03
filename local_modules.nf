//outputBarcode   = "${params.output}/barcode"
//outputDist      = "${params.output}/distribution"

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
	//publishDir params.output, mode: 'copy'
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
    //publishDir params.output, mode: 'copy'

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

process join_headers {

    container 'biocorecrg/doubletdetection:0.1'
    label (params.LABEL)
 	tag("${sampleID}")
	
	input:
	tuple val(sampleID), path("headers_*")

	output:
	tuple val(sampleID), path("header.sam")

	script:
	"""
		grep "@HD" headers_1 > header.sam 
		grep "@SQ" headers_1 >> header.sam 
		grep "@RG" headers_* --no-filename | huniq  >> header.sam 
 	"""


}


process make_distribution {
	//publishDir params.output, mode: 'copy'
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

process ectractAmpliconsFromPanelInfo {
        label (params.LABEL)
        tag("${panelinfo}")

        input:
        path(panelinfo)

        output:
        path("amplicons.csv")

        script:
        """
        cut -f 4 panel_info.bed | sort | uniq > amplicons.csv
        """
}

process cellfinder_scTAM {
        //publishDir params.output, mode: 'copy'
        container 'biocorecrg/tapestri_r:0.0.1'
        label (params.LABEL)
        tag("${sampleID}")

        input:
        tuple val(sampleID), path(input), path(amplifile)

        output:
        tuple val(sampleID), path("${sampleID}.barcode.cell.distribution.tsv")

        script:
        """
        cellfinder_scTAMseq.R -f ${input} -a ${amplifile} ${params.EXTRAPARS} -o ./${sampleID}.barcode.cell.distribution.tsv
        """

}

process cellfinder_MCL {
        //publishDir params.output, mode: 'copy'
        container 'biocorecrg/tapestri_r:0.0.1'
        label (params.LABEL)
        tag("${sampleID}")

        input:
        tuple val(sampleID), path(input), path(amplifile)

        output:
        tuple val(sampleID), path("${sampleID}.barcode.cell.distribution.tsv")

        script:
        """
        cellfinder_MCL_HC.R -f ${input} -a ${amplifile} ${params.EXTRAPARS} -o ./${sampleID}.barcode.cell.distribution.tsv
        """

}

process doubletDetection {

        container 'biocorecrg/doubletdetection:0.1'
        label (params.LABEL)
        tag("${sampleID}")

        input:
        tuple val(sampleID), path(barcode_distrib)

        output:
        tuple val(sampleID), path("${sampleID}_doublet_scores_DoubletDetection.csv")

        script:
        """
        DoubletDetection.py --input ${barcode_distrib} --output ${sampleID}_doublet_scores_DoubletDetection.csv
        """

}

process reads_per_cell {

        container 'biocorecrg/tapestri_r:0.0.2'
        label (params.LABEL)
        tag("${sampleID}")

        input:
        tuple val(sampleID), path(barcode_distrib), path(doublets), path(amplifile)

        output:
        tuple val(sampleID), path("${sampleID}_*")
 
        script:
        """
        plot_reads_per_cell.R --file ${barcode_distrib} --doublet ${doublets} \
                              --amplicon ${amplifile} --output ./ \
                              --name ${sampleID}
        """

}


process publishRes() {
        publishDir params.output, mode: 'copy'
        tag("${sampleID}")

        input:
        tuple val(sampleID), path(files)

        output:
        tuple val(sampleID), path("${sampleID}")


        script:
        """
		mkdir ${sampleID}
		mv ${files} ./${sampleID}
	        for i in ${sampleID}/*.txt; do cat \$i | gzip -c > ${sampleID}/`basename \$i`.gz; done
		rm ${sampleID}/*.txt      		
        """

}

// reshape the channel
def reshape(mychannel) {
	reshaped = mychannel.map{
		id = it[0].split("---")[0]
		[id, it[1]]
	}.groupTuple(sort: true)
	return(reshaped)
}

process protein_run {
	//publishDir params.output, mode: 'copy'
	container 'biocorecrg/tapestri:2.0.1'
    label (params.LABEL)
	tag("${sampleID}")
	
	input:
	tuple val(sampleID), path(reads), path(adapter3), path(adapter5), path(ab_barcodes)

	output:
	tuple val(sampleID), path("${sampleID}/results/tsv/${sampleID}-counts.tsv"), emit: out  

	script:
	
	def readsA = reads[0]
	def readsB = reads[1]
	
	"""
	tapestri protein run ${params.EXTRAPARS} \
		--ab-adapters-3 ${adapter3} \
		--ab-adapters-5 ${adapter5} \
		--ab-barcodes ${ab_barcodes} \
		--output-folder ./${sampleID} \
		--output-prefix ${sampleID} \
		--n-cores ${task.cpus} \
		--r1 ${readsA} \
		--r2 ${readsB}
	"""

}

process combine_prot_res {
	publishDir params.output, mode: 'copy'
	container 'biocorecrg/tapestri:2.0.1'
    label (params.LABEL)
	tag("${sampleID}")
	
	input:
	tuple val(sampleID), path(res)

	output:
	tuple val(sampleID), path("${sampleID}-counts.tsv"), emit: out  

	script:	
	"""
	combineRes.py -f "*.tsv" -o ${sampleID}-counts.tsv
	"""

}
