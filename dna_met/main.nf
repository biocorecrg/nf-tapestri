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

moduleFolder   = "${projectDir}/../BioNextflow/subworkflows"

include { final_message; notify_slack; colorCodes } from "${moduleFolder}/global_functions.nf"
def colors = colorCodes()


log.info """
${colors.yellow}${colors.bold}====================================================

╔╗ ┬┌─┐┌─┐┌─┐┬─┐┌─┐╔═╗╦═╗╔═╗  ╔╗╔╔═╗ ╔╦╗┌─┐┌─┐┌─┐┌─┐┌┬┐┬─┐┬  ╔╦╗╔╗╔╔═╗  ╔╦╗┌─┐┌┬┐
╠╩╗││ ││  │ │├┬┘├┤ ║  ╠╦╝║ ╦  ║║║╠╣───║ ├─┤├─┘├┤ └─┐ │ ├┬┘│   ║║║║║╠═╣  ║║║├┤  │ 
╚═╝┴└─┘└─┘└─┘┴└─└─┘╚═╝╩╚═╚═╝  ╝╚╝╚    ╩ ┴ ┴┴  └─┘└─┘ ┴ ┴└─┴  ═╩╝╝╚╝╩ ╩  ╩ ╩└─┘ ┴ 

====================================================

⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣶⣄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⣿⣸⣷⡄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢸⡏⣨⢿⣆⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣿⠉⢰⠿⣧⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢹⡗⠃⢀⡿⢷⡄⠀⠀⣠⣶⣦⣄⠀⢀⣤⣄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⢴⡆⠀⠀⠈⣿⡤⠞⠀⢸⢿⣦⡴⠟⠛⠓⢻⣷⢿⣵⢹⣆⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⢀⣽⣧⡀⠀⠀⢸⡆⢀⣴⡿⠋⠁⠀⠀⠀⣰⠟⡵⠋⢸⠈⣯⠛⢶⣤⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠙⠳⡾⠋⠁⠀⠀⣷⢋⡞⠁⢀⡴⠒⠋⣹⠟⡞⠁⠀⢸⠀⣿⠀⠀⠈⠻⣦⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⡇⠀⠀⠀⣠⡿⢿⡇⢰⡏⢀⣠⣼⡟⣸⠁⠀⠀⢸⠀⣿⠙⢦⡀⠀⠈⢷⣄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠉⠀⣠⡾⠋⠀⠈⢷⡈⢧⣸⣧⡘⠓⠇⠀⠀⠀⡞⠀⣿⠀⠀⠹⣄⠀⠈⡿⠳⣦⡀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣠⡾⠋⠀⠀⠀⠀⠀⠙⠓⠛⠋⠁⠀⠀⠀⠀⠸⠃⢸⡟⣧⠀⠀⠘⣆⣼⠇⠀⡾⣷⣀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⣠⡾⠋⠀⠀⠀⣀⣤⣤⣄⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⠀⢘⡷⠶⠖⠛⠁⠀⣴⣷⠏⠉⠻⣦⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⣠⣴⠟⠁⠀⠀⠀⢀⣾⣿⠉⣿⣿⣿⣦⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⢷⣀⣀⣤⣶⣿⡞⠉⠙⣆⠀⠹⣦⠀⠀⠀⠀
⠀⠀⢀⣴⠟⠉⡇⠀⠀⠀⠀⠀⠸⣿⣿⣿⣿⣤⣿⣽⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢹⠈⠈⡇⠀⡇⠀⠀⠘⡇⠀⢻⡆⠀⠀⠀
⠀⣰⢟⣑⡀⠸⡇⠀⠀⠀⠀⠀⠀⠻⣿⣿⣿⣻⡾⠋⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣼⠀⢀⡇⢠⡇⠀⠀⠀⢷⠀⠘⣇⠀⢀⡀
⢰⡏⠀⠿⠏⠀⢻⡄⠀⠀⠀⠀⠀⠀⠀⠉⠉⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⣀⡴⠃⢀⡞⣠⡟⢹⣄⠀⠀⢸⡀⠀⠻⠶⢾⡇
⢸⡇⠀⠀⠀⠀⠀⣿⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠛⠷⠶⣾⠷⠚⢹⡄⠀⠙⢦⣀⠀⠹⠦⠤⢤⣾⠃
⠀⢿⣄⠀⠀⠀⠀⢹⡄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢹⡄⠀⠀⠹⢦⡀⠀⠈⠛⠲⠶⡞⠻⣧⠀
⠀⠈⠻⢧⣄⣀⠀⣸⣇⣀⣀⣀⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⣴⠄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠻⣄⠀⠀⠀⠉⠓⢦⡀⠀⠀⣿⠀⢹⡇
⠀⠀⠀⠀⠈⠉⠉⠉⠁⠀⠈⠉⠙⠛⠶⢤⣤⣀⣀⣠⡤⠾⣻⠋⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢹⠳⣤⡀⠀⠀⠀⢹⡀⢀⡿⠀⢸⡇
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣸⠏⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣰⡟⠀⠀⢙⣷⡄⠀⢸⢁⣾⡇⠀⢸⡇
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢰⡟⠀⠀⠀⠀⠀⠀⠀⠀⣤⣤⠴⠞⠋⠀⢀⣠⠞⢀⣽⣠⡾⢛⡿⠀⢀⡿⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣿⠁⠀⠀⠀⠀⠀⠀⠀⠀⠈⠛⠿⠭⣿⣿⣭⣤⠴⠛⠙⠉⠀⡸⠁⢀⣾⠃⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢸⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⡟⠀⢀⡿⠀⠀⠀⠀⢰⡇⠀⡾⠁⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣾⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢸⠁⠀⢸⠀⠀⠀⠀⠀⠀⢳⣸⣇⣀⣀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢠⡿⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢸⡄⠀⠘⣆⠀⠀⠀⠀⠀⠀⠀⠉⣹⡇⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢸⣇⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⣷⠀⠀⠙⢦⣀⡀⠀⠀⣀⣠⣾⡟⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠉⠛⠶⠤⣤⣄⣀⣀⣀⣀⣀⣀⣀⣀⣤⡼⠷⣤⣀⡀⠈⠉⢉⣉⣭⠾⠋⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠉⠉⠉⠉⠉⠉⠉⠁⠀⠀⠀⠈⠉⠉⠛⠛⠋⠉⠁⠀⠀⠀⠀⠀
====================================================
${colors.reset}

${colors.bold}Input${colors.reset}
----------------------------------------------------${colors.reset}
${colors.green}reads${colors.reset}                        : ${params.reads}
${colors.green}genome${colors.reset}                       : ${params.genome}
${colors.green}panel${colors.reset}                        : ${params.panel}

${colors.bold}ampli$${colors.reset} can be '' in case of cellfinder_scTAMseq.R 
and it will take the 4th column of the panel file ${colors.reset}
----------------------------------------------------${colors.reset}
${colors.green}ampli${colors.reset}                        : ${params.ampli}

${colors.green}CpG_unmeth${colors.reset}                   : ${params.CpG_unmeth}

${colors.bold}cellfinder${colors.reset} can be ${colors.bold}cellfinder_scTAMseq.R${colors.reset} or ${colors.bold}cellfinder_MCL_HC.R${colors.reset}
----------------------------------------------------${colors.reset}
${colors.green}cellfinder${colors.reset}                  : ${params.cellfinder}

${colors.green}undigested${colors.reset}                   : ${params.undigested}
${colors.green}nreads${colors.reset}                       : ${params.nreads}
${colors.green}vdj${colors.reset}                          : ${params.vdj}
${colors.green}coord_f${colors.reset}                      : ${params.coord_f}
${colors.green}ref_trust${colors.reset}                    : ${params.ref_trust}
${colors.green}email${colors.reset}                        : ${params.email}
${colors.green}hook${colors.reset}                         : ${params.hook}
${colors.green}output (output folder)${colors.reset}       : ${params.output}
"""

if (params.resume) exit 1, "Are you making the classical --resume typo? Be careful!!!! ;)"

// INIZIALIZE MULTIQC REPORT
config_report = file("${projectDir}/config.yaml")
if( !config_report.exists() ) exit 1, "Missing config.yaml file!"
logo = file("${projectDir}/logo_dnamet.jpg")

localModule    = "${projectDir}/../local_modules.nf"
tools_opt      = "${projectDir}/tools_opt.tsv"


outputfolder     = "${params.output}"

include { getParameters } from "${localModule}"
progPars = getParameters(tools_opt)

include { barcode_extract } from "${localModule}" addParams(EXTRAPARS:progPars["tapestri-barcode-extract"], LABEL: "big_cpus")
include { make_target_counts } from "${localModule}" addParams(LABEL: "big_cpus")
include { make_minus_5panel } from "${localModule}"
include { make_distribution } from "${localModule}" addParams(LABEL: "big_mem_time")
include { doubletDetection } from "${localModule}" addParams(LABEL: "ignore")
include { ectractAmpliconsFromPanelInfo } from "${localModule}"
include { publishRes; reshape; reshapePairs; join_headers } from "${localModule}"
include { reformatPairA } from "${localModule}"addParams(LABEL: "big_mem")
include { reformatPairB } from "${localModule}" addParams(LABEL: "big_cpus")
include { extract_white_barcodes } from "${localModule}"

include { cellfinder_MCL } from "${localModule}" addParams( EXTRAPARS:progPars["cellfinder_MCL_HC.R"], LABEL: "big_mem" )
include { cellfinder_scTAM } from "${localModule}" addParams( EXTRAPARS:progPars["cellfinder_scTAMseq.R"], LABEL: "big_mem")
include { reads_per_cell } from "${localModule}" addParams ( LABEL: "big_mem" )
include { filter_cells } from "${localModule}" addParams ( LABEL: "big_mem_skip" )
include { cat_reads as cat_reads_A } from "${localModule}" addParams ( EXTRAPARS:"1" )
include { cat_reads as cat_reads_B } from "${localModule}" addParams ( EXTRAPARS:"2" )
include { cat_white_barcodes } from "${localModule}" 

include { trust_cluster } from "${localModule}" addParams ( EXTRAPARS:progPars["trust-cluster"] )
include { barcoderep_to_10X } from "${localModule}" 
include { join_metrics; cat_metrics } from "${localModule}" 

include { INDEX as BWA_INDEX } from "${moduleFolder}/alignment/bwa" addParams(STOREINDEX:params.index, CONTAINER:"biocorecrg/tapestri:2.0.1")
include { FASTQC_ID as FASTQC_QC } from "${moduleFolder}/qc/fastqc" 
include { REPORT as MULTIQC_REPORT } from "${moduleFolder}/reporting/multiqc" addParams(EXTRAPARS: "-c ${config_report.getName()}", OUTPUT: "${outputfolder}/report/")


include { CAT as SAMTOOLS_CAT } from "${moduleFolder}/misc/samtools" addParams(LABEL: "big_cpus")
include { SORT as SAMTOOLS_SORT } from "${moduleFolder}/misc/samtools" addParams(LABEL: "big_cpus")
include { STAT as SAMTOOLS_FLAGSTAT } from "${moduleFolder}/misc/samtools" addParams(LABEL: "big_cpus")
include { FASTQ_PAIRS_ONREADGROUP as SAMTOOLS_GET_UNMAPPED } from "${moduleFolder}/misc/samtools" addParams(GZIP:"YES", EXTRAPARS: " -f 1 -F 2 -T RG")


include { map_with_params as bwa_map } from "${moduleFolder}/alignment/bwa" addParams(CONTAINER:"biocorecrg/tapestri:2.0.1", EXTRAPARS:progPars["bwa"], LABEL: "big_cpus")
include { BVIEW as SAMTOOLS_VIEW } from "${moduleFolder}/misc/samtools" addParams(CONTAINER:"biocorecrg/tapestri:2.0.1", EXTRAPARS:progPars["samtools"], LABEL: "big_cpus")
include { BVIEW as SAMTOOLS_GET_HEADER } from "${moduleFolder}/misc/samtools" addParams(CONTAINER:"biocorecrg/tapestri:2.0.1", EXTRAPARS: " -H ", LABEL: "big_cpus")
include { FIXMATE_INFO as PICARD_FIXMATE } from "${moduleFolder}/misc/picard" addParams(CONTAINER:"biocorecrg/tapestri:2.0.1", EXTRAPARS:progPars["picard"], LABEL: "big_cpus")
include { INDEX as SAMTOOLS_INDEX } from "${moduleFolder}/misc/samtools" addParams(CONTAINER:"biocorecrg/tapestri:2.0.1", LABEL: "big_cpus")
include { RUNPE as RUNPE_TRUST4 } from "${moduleFolder}/immune_sys/trust4" addParams(EXTRAPARS:progPars["trust4"], LABEL: "big_time_cpus")


Channel.fromFilePairs( params.reads, flat: true) 
.ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
.set { ori_read_files }    

Channel.fromPath( params.panel, checkIfExists: true )
.ifEmpty { error "Cannot find any panel file ${params.panel}" }
.set { panel_file }    

undigested_par = channel.value(params.undigested)

Channel.from( config_report, logo ).set{data_for_reporting}


workflow {
	ori_read_files.map{
		[it[0], [it[1], it[2]] ]
	}.set{read_to_fastqc}
	fastqc_res = FASTQC_QC(read_to_fastqc).map{it[1]}.flatten()
	
    data_for_reporting = data_for_reporting.mix(fastqc_res)
    
    ori_read_files.splitFastq( by: params.nreads, pe: true, file: true, compress:true ).map{
        num = "${it[1]}".split("\\.")[-3]
        [ "${it[0]}---${num}", [ it[1], it[2] ] ]
    }.set{read_files}
   
    index = BWA_INDEX(params.genome)
    outbc = barcode_extract(read_files)

	// collectStats
	joint_metrics = join_metrics(reshape(outbc.barcode_metrics))
	
	all_joint_metrics = cat_metrics(joint_metrics.map{it[1]}.collect())

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
        
    data_for_reporting = data_for_reporting.mix(all_joint_metrics)
    
    //aln_stats = SAMTOOLS_FLAGSTAT(sorted_jaln)

    make_target_counts(index_reads)
    distros = make_distribution(index_reads.combine(m5panel))
    distrib = channel.empty()
    switch(params.cellfinder) {
        case "cellfinder_MCL_HC.R":
            Channel.fromPath( params.ampli, checkIfExists: true )
            .ifEmpty { error "Cannot find any ampli file ${params.ampli}" }.set { ampli_file }
            distrib = cellfinder_MCL(distros.dist_read1.combine(ampli_file))
            break;
        case "cellfinder_scTAMseq.R":
            if ( params.ampli == "") {
                ampli_file = ectractAmpliconsFromPanelInfo(panel_file) 
            } else {
                Channel.fromPath( params.ampli, checkIfExists: true )
                .ifEmpty { error "Cannot find any ampli file ${params.ampli}" }.set { ampli_file }      
            }
            distrib = cellfinder_scTAM(distros.dist_read1.combine(ampli_file))
            break;
    }
    
    doublets = doubletDetection(distrib)
    data_for_plot = distrib.join(doublets).combine(ampli_file)
        
    read_plots = reads_per_cell(data_for_plot)
    filtered_pdfs = channel.empty() 
    filtered_cells = channel.empty()
	trust_res = channel.empty()
	clus_res = channel.empty()
	t10x_res = channel.empty()

    if (params.CpG_unmeth != "") {

        Channel.fromPath ( params.CpG_unmeth, checkIfExists: true )
        .ifEmpty { error "Cannot find any CpG_unmeth file ${params.CpG_unmeth}" }
        .set { CpG_unmeth }

		data_for_filtering = distrib.join(doublets).combine(ampli_file).combine(CpG_unmeth).combine(undigested_par)
	
		filtered_res = filter_cells(data_for_filtering)  
		filtered_cells = filtered_res.filtered
		filtered_pdfs = filtered_res.pdfs
		good_bcs = filtered_res.good_bcs
		outmap.map{
			def sam_raw_ids = it[0].split("---") 
			[ sam_raw_ids[0], it[0], it[1] ] 
		}.combine(good_bcs, by: 0).map{
			[ it[1], it[2], it[3]]
		}.set{data_to_unmap}
		
		if (params.vdj == "YES") {

		// OLD FROM HERE
		/*
			bc_read_A = outbc.barcode_reads.map{
				[it[0], it[1]]
			}
	
			bc_read_B = outbc.barcode_reads.map{
				[it[0], it[2]]
			}    
		
		*/ 
		 // NEW HERE
           
            unmapped_fastqs = SAMTOOLS_GET_UNMAPPED(data_to_unmap)

			bc_read_A = unmapped_fastqs.map{
				[it[0], it[1]]
			}
	

			bc_read_B = unmapped_fastqs.map{
				[it[0], it[2]]
			}          
			
//		    ref_pairB = reformatPairB(bc_read_B)
		
//		    white_BCs = extract_white_barcodes(ref_pairA)

	 // COMMON HERE
	 	 	
		    ref_pairA = reformatPairA(bc_read_A)
		    ref_pairB = reformatPairB(bc_read_B)
		
		    white_BCs = extract_white_barcodes(bc_read_A)
 
 		    reformat_A = cat_reads_A(reshape(ref_pairA))
		    reformat_B = cat_reads_B(reshape(ref_pairB))
		    	  
		    white_BC = cat_white_barcodes(reshape(white_BCs))
		
		    ref_trust = channel.fromPath(params.ref_trust, checkIfExists: true)
		    coord_f = channel.fromPath(params.coord_f, checkIfExists: true)
		
		    reformat_reads = reformat_A.join(reformat_B)		    
            trust_res = RUNPE_TRUST4(ref_trust, coord_f, white_BC, reformat_reads)
            clus_res = trust_cluster(trust_res)
            t10x_res = barcoderep_to_10X(trust_res)
            
	    } 
    }
 
    
    MULTIQC_REPORT(data_for_reporting.collect())
    
    // collect results and publish using a custom output
    all_res = filtered_pdfs.transpose().mix(trust_res).mix(t10x_res.transpose()).mix(clus_res).mix(filtered_cells).mix(sorted_jaln).mix(indexes).mix(doublets).mix(distros.dist_read1).mix(distrib).mix(read_plots.transpose()).groupTuple().map{
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


