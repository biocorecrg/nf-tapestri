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

╔╗ ┬┌─┐┌─┐┌─┐┬─┐┌─┐╔═╗╦═╗╔═╗  ╔╦╗┌─┐┌─┐┌─┐┌─┐┌┬┐┬─┐┬  ╔═╗┬─┐┌─┐┌┬┐┌─┐┬┌┐┌┌─┐
╠╩╗││ ││  │ │├┬┘├┤ ║  ╠╦╝║ ╦   ║ ├─┤├─┘├┤ └─┐ │ ├┬┘│  ╠═╝├┬┘│ │ │ ├┤ ││││└─┐
╚═╝┴└─┘└─┘└─┘┴└─└─┘╚═╝╩╚═╚═╝   ╩ ┴ ┴┴  └─┘└─┘ ┴ ┴└─┴  ╩  ┴└─└─┘ ┴ └─┘┴┘└┘└─┘
                                                                                       
====================================================
reads                       : ${params.reads}
adapt3                      : ${params.adapt3}
adapt5                      : ${params.adapt5}
ab_barcode                  : ${params.ab_barcode}
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
include { protein_run } from "${localModule}" addParams(EXTRAPARS:progPars["tapestri-protein-run"], LABEL: "big_cpus")
include { combine_prot_res } from "${localModule}" 
include { publishRes; reshape } from "${localModule}"


Channel.fromFilePairs( params.reads, flat: true, checkIfExists: true) 
.ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
.set { ori_read_files }    

Channel.fromPath(params.adapt3, checkIfExists: true).set{adapt3}
Channel.fromPath(params.adapt5, checkIfExists: true).set{adapt5}
Channel.fromPath(params.ab_barcode, checkIfExists: true).set{ab_barcode}

workflow {
	ori_read_files.splitFastq( by: params.nreads, pe: true, file: true, compress:true ).map{
		num = "${it[1]}".split("\\.")[-3]
		[ "${it[0]}---${num}", [ it[1], it[2] ] ]
	}.set{read_files}

	prot_data = read_files.combine(adapt3).combine(adapt5).combine(ab_barcode)
    out_prot = protein_run(prot_data)
    comb_res = reshape(out_prot)
    combine_prot_res(comb_res)
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


