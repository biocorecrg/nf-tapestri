#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Copyright (c) 2023, Centre for Genomic Regulation (CRG) and the authors.

/*
===========================================================
 Nf-tapestri pipeline by Bioinformatics Core @ CRG

 @authors
 Luca Cozzuto <lucacozzuto@gmail.com>
=========================================================== 
*/

version = '0.1'

/*
*/

params.help            = false
params.resume          = false


log.info """
BIOCORE@CRG Tapestri - N F  ~  version ${version}

╔╗ ┬┌─┐┌─┐┌─┐┬─┐┌─┐╔═╗╦═╗╔═╗  ╔╗╔╔═╗ ╔╦╗┌─┐┌─┐┌─┐┌─┐┌┬┐┬─┐┬
╠╩╗││ ││  │ │├┬┘├┤ ║  ╠╦╝║ ╦  ║║║╠╣───║ ├─┤├─┘├┤ └─┐ │ ├┬┘│
╚═╝┴└─┘└─┘└─┘┴└─└─┘╚═╝╩╚═╚═╝  ╝╚╝╚    ╩ ┴ ┴┴  └─┘└─┘ ┴ ┴└─┴
                                                                                       
====================================================
reads                       : ${params.reads}
genome                      : ${params.genome}
email                       : ${params.email}
output (output folder)      : ${params.output}
"""

if (params.resume) exit 1, "Are you making the classical --resume typo? Be careful!!!! ;)"

moduleFolder   = "${projectDir}/BioNextflow/subworkflows"
multiQC_conf   = Channel.fromPath("${projectDir}/config.yaml")
localModule    = "${projectDir}/local_modules.nf"
tools_opt      = "${projectDir}/tools_opt.tsv"
//adapters       = "${projectDir}/adapters.txt"

//if( !genome_file.exists() ) exit 1, "Missing genome file: ${genome_file}"

outputfolder    = "${params.output}"

include { getParameters } from "${localModule}"
progPars = getParameters(tools_opt)

include { final_message; notify_slack } from "${moduleFolder}/global_functions.nf"
include { barcode_extract } from "${localModule}" addParams(EXTRAPARS:progPars["tapestri-barcode-extract"], LABEL: "big_cpus")
include { INDEX as BWA_INDEX } from "${moduleFolder}/alignment/bwa" addParams(CONTAINER:"biocorecrg/tapestri:2.0.1")
include { map_with_params as bwa_map } from "${moduleFolder}/alignment/bwa" addParams(CONTAINER:"biocorecrg/tapestri:2.0.1", EXTRAPARS:progPars["bwa"], LABEL: "big_cpus")
include { BVIEW as SAMTOOLS_VIEW } from "${moduleFolder}/misc/samtools" addParams(CONTAINER:"biocorecrg/tapestri:2.0.1", EXTRAPARS:progPars["samtools"], LABEL: "big_cpus")

Channel.fromFilePairs( params.reads, size: 2 )
.ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
.set { read_files }    



workflow {
	outbc = barcode_extract(read_files)
	index = BWA_INDEX(params.genome)
	outbc.barcode_reads.join(outbc.sam_header).map{
		[ it[0], [ it[1], it[2] ], "-H ${it[3]}" ]
	}.set{ data_for_alignment }
	outmap = bwa_map(data_for_alignment, index)
	outmap.view()
	SAMTOOLS_VIEW(outmap)

}



workflow.onComplete {

    def text = final_message("NF-Tapestri")
    println text
    if (params.hook != "") {
    	notify_slack(text, params.hook)
	}
}

/*
* Mail notification
*/

if (params.email == "yourmail@yourdomain" || params.email == "") { 
    log.info 'Skipping the email\n'
}
else {
    log.info "Sending the email to ${params.email}\n"

    workflow.onComplete {
 	   def msg = final_message("NF-Tapestri")	
        sendMail(to: params.email, subject: "NF-Tapestri execution", body: msg, attach: "${outputMultiQC}/multiqc_report.html")
    }
}


