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
panel                       : ${params.panel}
ampli                       : ${params.ampli}
cellfinder                  : ${params.cellfinder}
email                       : ${params.email}
output (output folder)      : ${params.output}
"""

if (params.resume) exit 1, "Are you making the classical --resume typo? Be careful!!!! ;)"

moduleFolder   = "${projectDir}/BioNextflow/subworkflows"
multiQC_conf   = Channel.fromPath("${projectDir}/config.yaml")
localModule    = "${projectDir}/local_modules.nf"
tools_opt      = "${projectDir}/tools_opt.tsv"
//adapters       = "${projectDir}/adapters.txt"


outputfolder     = "${params.output}"
outputAlignments = "${params.output}/bam"

include { getParameters } from "${localModule}"
progPars = getParameters(tools_opt)

include { final_message; notify_slack } from "${moduleFolder}/global_functions.nf"
include { barcode_extract } from "${localModule}" addParams(EXTRAPARS:progPars["tapestri-barcode-extract"], LABEL: "big_cpus")
include { make_target_counts } from "${localModule}" addParams(LABEL: "big_cpus")
include { make_minus_5panel } from "${localModule}"
include { make_distribution } from "${localModule}"
include { ectractAmpliconsFromPanelInfo } from "${localModule}"
include { publishRes } from "${localModule}"

include { cellfinder_MCL } from "${localModule}" addParams( EXTRAPARS:progPars["cellfinder_MCL_HC.R"] )
include { cellfinder_scTAM } from "${localModule}" addParams( EXTRAPARS:progPars["cellfinder_scTAMseq.R"] )

include { INDEX as BWA_INDEX } from "${moduleFolder}/alignment/bwa" addParams(CONTAINER:"biocorecrg/tapestri:2.0.1")
include { map_with_params as bwa_map } from "${moduleFolder}/alignment/bwa" addParams(CONTAINER:"biocorecrg/tapestri:2.0.1", EXTRAPARS:progPars["bwa"], LABEL: "big_cpus")
include { BVIEW as SAMTOOLS_VIEW } from "${moduleFolder}/misc/samtools" addParams(CONTAINER:"biocorecrg/tapestri:2.0.1", EXTRAPARS:progPars["samtools"], LABEL: "big_cpus")
include { FIXMATE_INFO as PICARD_FIXMATE } from "${moduleFolder}/misc/picard" addParams(CONTAINER:"biocorecrg/tapestri:2.0.1", EXTRAPARS:progPars["picard"], LABEL: "big_cpus")
include { INDEX as SAMTOOLS_INDEX } from "${moduleFolder}/misc/samtools" addParams(CONTAINER:"biocorecrg/tapestri:2.0.1", LABEL: "big_cpus")

Channel.fromFilePairs( params.reads, size: 2 )
.ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
.set { read_files }    

Channel.fromPath( params.panel )
.ifEmpty { error "Cannot find any panel file ${params.panel}" }
.set { panel_file }    


workflow {
    outbc = barcode_extract(read_files)
    index = BWA_INDEX(params.genome)
    outbc.barcode_reads.join(outbc.sam_header).map{
        [ it[0], [ it[1], it[2] ], "-H ${it[3]}" ]
    }.set{ data_for_alignment }
    outmap = bwa_map(data_for_alignment, index)
    bams = SAMTOOLS_VIEW(outmap)
    fixed = PICARD_FIXMATE(bams)
    indexes = SAMTOOLS_INDEX(fixed)
    m5panel = make_minus_5panel(panel_file)
    index_reads = fixed.join(indexes) 
    make_target_counts(index_reads)
    distros = make_distribution(index_reads.combine(m5panel))

    distrib = channel.empty()
    if ( params.ampli == "") {
        ampli_file = ectractAmpliconsFromPanelInfo(panel_file) 
    } else {
         Channel.fromPath( params.ampli )
         .ifEmpty { error "Cannot find any ampli file ${params.ampli}" }.set { ampli_file }
    }

    switch(params.cellfinder) {
        case "cellfinder_MCL_HC.R":
            distrib = cellfinder_MCL(distros.dist_read1.combine(ampli_file))
            break;
        case "cellfinder_scTAMseq.R":
            distrib = cellfinder_scTAM(distros.dist_read1.combine(ampli_file))
            break;
    }
    // collect results and publish using a custom output
    all_res = fixed.mix(indexes).mix(distros.dist_read1).mix(distrib).groupTuple().map{
		def name = it[0].split('_')[0]
    	[name, it[1]]
    }
    publishRes(all_res)
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


