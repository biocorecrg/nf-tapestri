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
nreads                      : ${params.nreads}
email                       : ${params.email}
output (output folder)      : ${params.output}
"""

if (params.resume) exit 1, "Are you making the classical --resume typo? Be careful!!!! ;)"

moduleFolder   = "${projectDir}/../BioNextflow/subworkflows"
multiQC_conf   = Channel.fromPath("${projectDir}/config.yaml")
localModule    = "${projectDir}/../local_modules.nf"
tools_opt      = "${projectDir}/tools_opt.tsv"


outputfolder     = "${params.output}"

include { getParameters } from "${localModule}"
progPars = getParameters(tools_opt)

include { final_message; notify_slack } from "${moduleFolder}/global_functions.nf"
include { barcode_extract } from "${localModule}" addParams(EXTRAPARS:progPars["tapestri-barcode-extract"], LABEL: "big_cpus")
include { make_target_counts } from "${localModule}" addParams(LABEL: "big_cpus")
include { make_minus_5panel } from "${localModule}"
include { make_distribution } from "${localModule}" addParams(LABEL: "big_mem_time")
include { doubletDetection } from "${localModule}" addParams(LABEL: "ignore")
include { ectractAmpliconsFromPanelInfo } from "${localModule}"
include { publishRes; reshape; join_headers } from "${localModule}"

include { cellfinder_MCL } from "${localModule}" addParams( EXTRAPARS:progPars["cellfinder_MCL_HC.R"], LABEL: "big_mem" )
include { cellfinder_scTAM } from "${localModule}" addParams( EXTRAPARS:progPars["cellfinder_scTAMseq.R"], LABEL: "big_mem")
include { reads_per_cell } from "${localModule}" addParams ( LABEL: "big_mem" )

include { INDEX as BWA_INDEX } from "${moduleFolder}/alignment/bwa" addParams(STOREINDEX:"${outputfolder}/index", CONTAINER:"biocorecrg/tapestri:2.0.1")
include { FASTQCP as FASTQC_QC } from "${moduleFolder}/qc/fastqc" 
include { REPORT as MULTIQC_REPORT  } from "${moduleFolder}/reporting/multiqc" addParams(OUTPUT: "${outputfolder}/report/")

include { CAT as SAMTOOLS_CAT } from "${moduleFolder}/misc/samtools" addParams(LABEL: "big_cpus")
include { SORT as SAMTOOLS_SORT } from "${moduleFolder}/misc/samtools" addParams(LABEL: "big_cpus")
include { STAT as SAMTOOLS_FLAGSTAT } from "${moduleFolder}/misc/samtools" addParams(LABEL: "big_cpus")


include { map_with_params as bwa_map } from "${moduleFolder}/alignment/bwa" addParams(CONTAINER:"biocorecrg/tapestri:2.0.1", EXTRAPARS:progPars["bwa"], LABEL: "big_cpus")
include { BVIEW as SAMTOOLS_VIEW } from "${moduleFolder}/misc/samtools" addParams(CONTAINER:"biocorecrg/tapestri:2.0.1", EXTRAPARS:progPars["samtools"], LABEL: "big_cpus")
include { BVIEW as SAMTOOLS_GET_HEADER } from "${moduleFolder}/misc/samtools" addParams(CONTAINER:"biocorecrg/tapestri:2.0.1", EXTRAPARS: " -H ", LABEL: "big_cpus")
include { FIXMATE_INFO as PICARD_FIXMATE } from "${moduleFolder}/misc/picard" addParams(CONTAINER:"biocorecrg/tapestri:2.0.1", EXTRAPARS:progPars["picard"], LABEL: "big_cpus")
include { INDEX as SAMTOOLS_INDEX } from "${moduleFolder}/misc/samtools" addParams(CONTAINER:"biocorecrg/tapestri:2.0.1", LABEL: "big_cpus")

Channel.fromFilePairs( params.reads, flat: true) 
.ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
.set { ori_read_files }    

Channel.fromPath( params.panel )
.ifEmpty { error "Cannot find any panel file ${params.panel}" }
.set { panel_file }    


workflow {
    data_for_reporting = FASTQC_QC(ori_read_files)
    
    ori_read_files.splitFastq( by: params.nreads, pe: true, file: true, compress:true ).map{
        num = "${it[1]}".split("\\.")[-3]
        [ "${it[0]}---${num}", [ it[1], it[2] ] ]
    }.set{read_files}
    
    index = BWA_INDEX(params.genome)
    outbc = barcode_extract(read_files)

    outbc.barcode_reads.join(outbc.sam_header).map{
        [ it[0], [ it[1], it[2] ], "-H ${it[3]}" ]
    }.set{ data_for_alignment } 
    outmap = bwa_map(data_for_alignment, index)
    bams = SAMTOOLS_VIEW(outmap)
    fixed = PICARD_FIXMATE(bams)
    fixed_headers = SAMTOOLS_GET_HEADER(fixed)
    sam_header = join_headers(reshape(fixed_headers))
    data_2_cat = reshape(fixed).join(sam_header)

    jaln_reads = SAMTOOLS_CAT(data_2_cat)
    sorted_jaln = SAMTOOLS_SORT(jaln_reads) 
    indexes = SAMTOOLS_INDEX(sorted_jaln)

    m5panel = make_minus_5panel(panel_file)
    index_reads = sorted_jaln.join(indexes) 
    
    //aln_stats = SAMTOOLS_FLAGSTAT(sorted_jaln)

    make_target_counts(index_reads)
    distros = make_distribution(index_reads.combine(m5panel))
    distrib = channel.empty()
    switch(params.cellfinder) {
        case "cellfinder_MCL_HC.R":
            Channel.fromPath( params.ampli )
            .ifEmpty { error "Cannot find any ampli file ${params.ampli}" }.set { ampli_file }
            distrib = cellfinder_MCL(distros.dist_read1.combine(ampli_file))
            break;
        case "cellfinder_scTAMseq.R":
            if ( params.ampli == "") {
                ampli_file = ectractAmpliconsFromPanelInfo(panel_file) 
            } else {
                Channel.fromPath( params.ampli )
                .ifEmpty { error "Cannot find any ampli file ${params.ampli}" }.set { ampli_file }      
            }
            distrib = cellfinder_scTAM(distros.dist_read1.combine(ampli_file))
            break;
    }
    
    doublets = doubletDetection(distrib)
    data_for_plot = distrib.join(doublets).combine(ampli_file)
        
    read_plots = reads_per_cell(data_for_plot)
    
    MULTIQC_REPORT(data_for_reporting.collect())
    
    // collect results and publish using a custom output
    all_res = sorted_jaln.mix(indexes).mix(doublets).mix(distros.dist_read1).mix(distrib).mix(read_plots.transpose()).groupTuple().map{
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


