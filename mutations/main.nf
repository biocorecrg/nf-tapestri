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
genome                      : ${params.genome}
gatk_config                 : ${params.gatk_config}
panel                       : ${params.panel}
amplicons                   : ${params.amplicons}
email                       : ${params.email}
output (output folder)      : ${params.output}
"""

if (params.resume) exit 1, "Are you making the classical --resume typo? Be careful!!!! ;)"

moduleFolder   = "${projectDir}/../BioNextflow/subworkflows"
multiQC_conf   = Channel.fromPath("${projectDir}/config.yaml")
localModule    = "${projectDir}/../local_modules.nf"
tools_opt      = "${projectDir}/tools_opt.tsv"
outputfolder   = "${params.output}"

include { getParameters } from "${localModule}"
progPars = getParameters(tools_opt)


include { final_message; notify_slack } from "${moduleFolder}/global_functions.nf"
include { cell_bam } from localModule addParams(LABEL: 'big_cpus') 
//include { genotyping_split  } from localModule addParams(LABEL: 'big_cpus') 
include { INDEX as SAMTOOLS_INDEX } from "${moduleFolder}/misc/samtools" addParams(LABEL: "big_cpus")
include { genotyping_split_and_haplotype } from localModule addParams(LABEL: 'big_mem_time_cpus') 

include { changeGATKconfig } from "${localModule}"
include { PossiblyUnzipGenome } from "${moduleFolder}/misc/misc.nf"
include { FAIDX as SAMTOOLS_FAIDX } from "${moduleFolder}/misc/samtools" addParams(CONTAINER:"biocorecrg/tapestri:2.0.1")
include { DICT as SAMTOOLS_DICT } from "${moduleFolder}/misc/samtools" addParams(CONTAINER:"biocorecrg/tapestri:2.0.1")

include { genotyping_import } from localModule addParams(LABEL: 'big_mem') 
include { genotyping } from localModule addParams(EXTRAPARS:progPars["genotyping"], LABEL: "big_mem_time")
include { create_loom } from localModule addParams(LABEL: "big_mem_time")


inputFiles = channel.fromPath(params.input, checkIfExists: true)
genome = channel.fromPath(params.genome, checkIfExists: true)
gatk_config = channel.fromPath(params.gatk_config, checkIfExists: true)
panel = channel.fromPath(params.panel, checkIfExists: true)
amplicons = channel.fromPath(params.amplicons, checkIfExists: true)

inputFiles.branch {
        bam: it =~ /\.bam$/
        	return [it.getParent().getParent().getName(), it]
        bai: it =~ /\.bai$/
            return [it.getParent().getParent().getName(), it]
        cells: it =~ /filtered.tsv$/
            return [it.getParent().getParent().getName(), it]
        metadata: it =~ /metadata.csv$/
            return [it.getParent().getParent().getName(), it]
    }.set { input_data }


data_for_cell_bam = input_data.cells.join(input_data.bam).join(input_data.bai)

data_for_cell_bam.view()


workflow {

    checked_genome = PossiblyUnzipGenome(params.genome)
    faidx = SAMTOOLS_FAIDX(checked_genome)
    dict = SAMTOOLS_DICT(checked_genome)
    bam_cells = cell_bam(data_for_cell_bam)
    bai_cells = SAMTOOLS_INDEX(bam_cells)
    bams_indx = bam_cells.join(bai_cells)
    vcf_dir = genotyping_split_and_haplotype(bams_indx.combine(gatk_config).combine(checked_genome).combine(panel).combine(amplicons).combine(faidx).combine(dict))
    geno_out = genotyping_import(vcf_dir.combine(panel))
	vcf = genotyping(geno_out.combine(checked_genome).combine(faidx).combine(dict).combine(panel))
    vcf_data = vcf.join(input_data.metadata).combine(panel)
    vcf_data.view()
//	create_loom(vcf_data)
	
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


