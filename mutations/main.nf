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

╔╗ ┬┌─┐┌─┐┌─┐┬─┐┌─┐╔═╗╦═╗╔═╗  ╔╗╔╔═╗ ╔╦╗┌─┐┌─┐┌─┐┌─┐┌┬┐┬─┐┬  ╔╦╗┬ ┬┌┬┐┌─┐┌┬┐┬┌─┐┌┐┌┌─┐
╠╩╗││ ││  │ │├┬┘├┤ ║  ╠╦╝║ ╦  ║║║╠╣───║ ├─┤├─┘├┤ └─┐ │ ├┬┘│  ║║║│ │ │ ├─┤ │ ││ ││││└─┐
╚═╝┴└─┘└─┘└─┘┴└─└─┘╚═╝╩╚═╚═╝  ╝╚╝╚    ╩ ┴ ┴┴  └─┘└─┘ ┴ ┴└─┴  ╩ ╩└─┘ ┴ ┴ ┴ ┴ ┴└─┘┘└┘└─┘
                                                                                       
====================================================
input                       : ${params.input}
gatk_config                 : ${params.gatk_config}
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
include { cell_bam } from localModule addParams(LABEL: 'big_cpus') 
include { genotyping_split } from localModule addParams(LABEL: 'big_cpus') 
include { INDEX as SAMTOOLS_INDEX } from "${moduleFolder}/misc/samtools" addParams(LABEL: "big_cpus")

inputFiles = channel.fromPath(params.input, checkIfExists: true)
gatk_config = channel.fromPath(params.gatk_config, checkIfExists: true)

inputFiles.branch {
        bam: it =~ /\.bam$/
        	return [it.getParent().getParent().getName(), it]
        bai: it =~ /\.bai$/
            return [it.getParent().getParent().getName(), it]
        cells: it =~ /filtered.tsv$/
            return [it.getParent().getParent().getName(), it]
    }.set { input_data }

data_for_cell_bam = input_data.cells.join(input_data.bam).join(input_data.bai).view()


workflow {
	bam_cells = cell_bam(data_for_cell_bam)
	bai_cells = SAMTOOLS_INDEX(bam_cells)
	sam_dir = genotyping_split(bam_cells.join(bai_cells))
	genotyping_haplotype(sam_dir.join(gatk_config))

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


