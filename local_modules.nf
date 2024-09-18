params.EXTRAPARS = ""
params.LABEL = ""
params.OUTPUT = ""

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

process filter_cells {

        container 'biocorecrg/tapestri_r:0.0.3'
        label (params.LABEL)
        tag("${sampleID}")

        input:
        tuple val(sampleID), path(barcode_distrib), path(doublets), path(amplifile), path(CpG_unmeth) 

        output:
        tuple val(sampleID), path("${sampleID}.barcode.cell.distribution_filtered.tsv"), emit: filtered
        tuple val(sampleID), path("*.pdf"), emit: pdfs
 
        script:
        """
             create.object.QC.MB.pipeline.R ${sampleID} \
             ${barcode_distrib} \
             ${doublets} \
             ./ \
             ${amplifile} \
             ${CpG_unmeth} 
        """

}

process reformatPairA {

    container 'quay.io/biocontainers/biopython:1.68--py35_0'
        label (params.LABEL)
        tag("${sampleID}")

        input:
        tuple val(sampleID), path(fastqA)

        output:
        tuple val(sampleID), path("*_reformat.fastq.gz")
         
        script:
        def basefile = "${fastqA}".replaceAll('.with.barcodes_reformat.fastq.gz', '') 

        """
        reformat.py --listOfFastqs ${fastqA} --sampleName ${basefile} --outputFolder ./ 
        """

}

process reformatPairB {

        container 'quay.io/biocontainers/pigz:2.8'
        label (params.LABEL)
        tag("${sampleID}")

        input:
        tuple val(sampleID), path(fastqA)

        output:
        tuple val(sampleID), path("*_reformat.fastq.gz")
         
        script:
        def basefile = "${fastqA}".replaceAll('.with.barcodes_reformat.fastq.gz', '') 

        """
        pigz -dkc ${fastqA} | awk 'NR%4==1 {print \$1,\$3,"2:N:0:1"} NR%4!=1 {print \$0}' | pigz -p ${task.cpus} -c > ./${basefile}.with.barcodes_reformat.fastq.gz 
        """

}

process cat_reads {

        container 'quay.io/biocontainers/pigz:2.8'
        label (params.LABEL)
        tag("${sampleID}")

        input:
        tuple val(sampleID), path(fastqs)

        output:
        tuple val(sampleID), path("*.with.barcodes_reformat.fastq.gz")
         
        script:

        """
        cat ${fastqs} >> ${sampleID}_read${params.EXTRAPARS}.with.barcodes_reformat.fastq.gz
        """

}

process cat_white_barcodes {

        container 'quay.io/biocontainers/pigz:2.8'
        label (params.LABEL)
        tag("${sampleID}")

        input:
        tuple val(sampleID), path(whitebc)

        output:
        tuple val(sampleID), path("*_read1_MCL_barcodes.txt.gz")
         
        script:

        """
        cat ${whitebc} >> ${sampleID}_read1_MCL_barcodes.txt.gz
        """

}

process extract_white_barcodes {

        container 'quay.io/biocontainers/pigz:2.8'
        label (params.LABEL)
        tag("${sampleID}")

        input:
        tuple val(sampleID), path(fastq)

        output:
        tuple val(sampleID), path("*_MCL_barcodes.txt.gz")
         
        script:
        def basefile = "${fastq}".replaceAll('.with.barcodes_reformat.fastq.gz', '') 

        """
        zcat ${fastq} | grep '@'| cut -d " " -f 2 | sed 's/RG:Z://g'| sed 's/-1//g' | pigz -c -p ${task.cpus} > ${basefile}_MCL_barcodes.txt.gz 
        """

}

process trust_cluster {

        container 'quay.io/biocontainers/biopython:1.68--py35_0'
        label (params.LABEL)
        tag("${sampleID}")

        input:
        tuple val(sampleID), path(trust4_output)

        output:
        tuple val(sampleID), path("${trust4_output}_cdr3_clones.out")
         
        script:

        """
        trust-cluster.py ${trust4_output}/${trust4_output}_read1_cdr3.out ${params.EXTRAPARS} > ${trust4_output}_cdr3_clones.out
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
		mkdir ${sampleID}/output
		mkdir ${sampleID}/intermediate
		mkdir ${sampleID}/qc_log
		mkdir ${sampleID}/plots
		mkdir ${sampleID}/trust4
		mv *.pdf ${sampleID}/plots
		mv *.bam* ./${sampleID}/intermediate
		mv *_stats.txt ./${sampleID}/qc_log
		mv *_metadata.csv ./${sampleID}/qc_log
		mv *.png ./${sampleID}/qc_log
		mv *.txt ./${sampleID}/intermediate
		mv *.tsv ./${sampleID}/output
		mv *.csv ./${sampleID}/output
		for i in ${sampleID}/intermediate/*.txt; do cat \$i | gzip -c > ${sampleID}/intermediate/`basename \$i`.gz; done
		rm ${sampleID}/intermediate/*.txt
		if compgen -G "${sampleID}_*" > /dev/null; then
			mv ${sampleID}_* ${sampleID}/trust4
		fi
        """

}

/*
*/

// reshape the channel
def reshape(mychannel) {
	reshaped = mychannel.map{
		id = it[0].split("---")[0]
		[id, it[1]]
	}.groupTuple(sort: true)
	return(reshaped)
}

// reshape the channel
def reshapePairs(mychannel) {
	pairA = mychannel.map{
		id = it[0].split("---")[0]
		[id, it[1]]
	}.groupTuple(sort: true)

	pairB = mychannel.map{
		id = it[0].split("---")[0]
		[id, it[2]]
	}.groupTuple(sort: true)

    reshaped = pairA.join(pairB)
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

process cell_bam {
//	publishDir params.OUTPUTDIR, mode: 'copy'
	container 'biocorecrg/tapestri:2.0.1'
    label (params.LABEL)
	tag("${sampleID}")
	
	input:
	tuple val(sampleID), path(input_cells), path(bam), path(bamid)

	output:
	tuple val(sampleID), path("${sampleID}.cells.bam"), emit: out  

	script:	
	"""
    awk '{if (NR > 1) print \$1}' ${input_cells} > cells.txt
    samtools view -@ ${task.cpus} -b -R cells.txt ${bam} > ${sampleID}.cells.bam
	"""

}



process genotyping_split_and_haplotype {

	container 'biocorecrg/tapestri:2.0.1'
    label (params.LABEL)
	tag("${sampleID}")
	
	input:
	tuple val(sampleID), path(bam), path(bamid), path(gatkconfig_path), path(genome), path("panel_info.bed"), path("panel_info.amplicons"), path(faidx), path(dict)

	output:
	tuple val(sampleID), path("${sampleID}_vcf"), emit: out  

	script:	
	"""
       split_bam_by_barcode.py -i ${bam} -o ${sampleID}_bam -n ${task.cpus}
       
      sed s@'REPLACE_GENOME_PATH'@${genome}@ ${gatkconfig_path} | \
      sed s@'REPLACE_PANEL_INFO_BED'@panel_info.bed@  | \
      sed --expression s@'REPLACE_MCL_PATH'@\$PWD@ > gatk_new_config.yaml
      
      gatk_haplotype_caller.py --input-dir ${sampleID}_bam --output-dir ${sampleID}_vcf --cores ${task.cpus} --config gatk_new_config.yaml --prefix gatk

	"""

}

process genotyping_haplotype {

	container 'biocorecrg/tapestri:2.0.1'
    label (params.LABEL)
	tag("${sampleID}")
	
	input:
	tuple val(sampleID), path(inputdir), path(gatkconfig_path), path(genome), path(panelinfo), path(faidx), path(dict)

	output:
	tuple val(sampleID), path("${sampleID}_vcf"), emit: out  

	script:	
    
	"""
      sed s@'REPLACE_GENOME_PATH'@${genome}@ ${gatkconfig_path} | \
      sed s@'REPLACE_PANEL_INFO_BED'@${panelinfo}@  | \
      sed --expression s@'REPLACE_MCL_PATH'@\$PWD@ > gatk_new_config.yaml
      ln -s ${panelinfo} ${panelinfo.baseName}.amplicons 
      gatk_haplotype_caller.py --input-dir ${inputdir} --output-dir ${sampleID}_vcf --cores ${task.cpus} --config gatk_new_config.yaml --prefix gatk

	"""

}

/*
*/
process genotyping_import {

	container 'biocorecrg/tapestri:2.0.1'
    label (params.LABEL)
	tag("${sampleID}")
	
	input:
	tuple val(sampleID), path(vcf_folder), path(panel)

	output:
	tuple val(sampleID), path("${sampleID}_gendb"), emit: out  

	script:
    
	"""
        call_string=""
        for f in ${vcf_folder}/*.vcf.gz; do
                call_string+=" --variant \$f"
        done ;
        gatk GenomicsDBImport \$call_string --genomicsdb-workspace-path ${sampleID}_gendb --intervals ${panel} --merge-input-intervals --batch-size 150  
	"""

}

/*
*/
process genotyping {

	container 'biocorecrg/tapestri:2.0.1'
    label (params.LABEL)
	tag("${sampleID}")
	
	input:
	tuple val(sampleID), path(db_folder), path(genome), path(fai), path(dict), path(panel)

	output:
	tuple val(sampleID), path("${sampleID}.combined.vcf.gz"), emit: out  

	script:
    
	"""
	cp /nfs/scratch01/bi/luca/tap_workaround/${sampleID}.combined.vcf .
	bgzip ${sampleID}.combined.vcf	

 	"""
//        export JAVA_OPTIONS="-Xmx${task.memory.toGiga()/2}g"
//        alias java="java -Xmx${task.memory.toGiga()/2}g" ;
//        gatk GenotypeGVCFs --java-options "-Xmx${task.memory.toGiga()/2}g -XX:ParallelGCThreads=${task.cpus}" -R ${genome} ${params.EXTRAPARS} -V gendb://${db_folder} -L ${panel} -O ${sampleID}.combined.vcf

}

/*
*/
process create_loom {

	publishDir params.OUTPUT, mode: 'copy'

	container 'biocorecrg/tapestri:2.0.1'
    label (params.LABEL)
	tag("${sampleID}")
	
	input:
	tuple val(sampleID), path(gzvcf), path(metadata), path(panel_bed) 

	output:
	tuple val(sampleID), path("${sampleID}.dna.loom")

	script:
	"""
	    mkdir tmp
        printf "genome:\n %s version: hg19" > ./tmp/hg19.yaml ;
        vcf_to_loom.py \
        ./tmp/hg19.yaml \
        ${panel_bed} \
        ${metadata} \
        <(zcat ${gzvcf}) \
        ${sampleID}.dna.loom \
        --logfile log.txt 
 	"""

}

process create_h5 {
	publishDir params.OUTPUT, mode: 'copy'

	container 'biocorecrg/tapestri:2.0.1'
    label (params.LABEL)
	tag("${sampleID}")
	
	input:
	tuple val(sampleID), path(bgvcf), path(cells), path(loom) 

	output:
	tuple val(sampleID), path("${sampleID}.dna.h5")

	script:
	"""
        printf "{\n \\"sample_name\\": \\"${sampleID}\\",\n \\"pipeline_version\\": \\"2.0.1\\",\n \\"genome_version\\": \\"hg19\\",\n \\"panel_name\\": \\"MCL\\",\n \\"n_amplicons\\": 910,\n \\"chemistry_version\\": \\"V2\\"\n}" > metadata.json ;
        tapestri h5 create dna \
        --vcf ${bgvcf} \
        --read-counts ${cells} \
        --metadata metadata.json \
        --output ${sampleID}.dna.h5 ;
 	"""


}



process changeGATKconfig {

	//publishDir params.output, mode: 'copy'
	container 'biocorecrg/tapestri:2.0.1'
	
	input:
	tuple path(gatkconfig_path), path(genome), val(panelinfo)

	output:
	path("gatk_new_config.yaml")

	script:
	def panel_path = file(panelinfo).getParent()
	
	"""
       sed s@'REPLACE_GENOME_PATH'@${genome}@ ${gatkconfig_path} | \
       sed s@'REPLACE_PANEL_INFO_BED'@${panelinfo}@ > gatk_new_config.yaml
	"""
}




