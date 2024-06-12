#!/usr/bin/env python

import argparse
import io
import logging
import threading
import time
import yaml

from missionbio.core.createMetaforLoom import create_meta_data
from missionbio.core.log_utils import configure_root_logger
from missionbio.core.loomifier import (TapData, LoomWriter)
from missionbio.core.vcf import MBVcf


logger = logging.getLogger(__name__)

REQUIRED_FIELDS_IN_VCF = [
    'CHROM', 'POS', 'GT', 'calldata/DP', 'AD', 'GQ', 'samples', 'INFO',
    'ID', 'REF', 'ALT', 'QUAL']


def thread_task(task, amp_bedfile, logger):

    global mb_vcf
    global data_dict
    global threadLock

    try:
        if task == 'NGT':
            mb_vcf.numerical_genotypebypos()
        elif task == 'DP':
            mb_vcf.readdepthbypos()
        elif task == 'GQ':
            mb_vcf.genotypequalitybypos()
        elif task == 'AD':
            mb_vcf.allelicdepthbypos()
        elif task == 'RO':
            mb_vcf.referencealleles()

        task_result = mb_vcf.filter_by_quality(task)
        task_result = mb_vcf.replace_na_with_zero(task_result)

        threadLock.acquire()
        data_dict[task] = task_result
        logger.info("{}: {}".format(task, data_dict[task]))
        threadLock.release()

    except Exception as e:
        msg = "Thread {} failed due to {}".format(
            threading.get_ident(), str(e))
        logger.error(msg)
        raise Exception(msg)


def prepare_snp_info(task_result, amp_bedfile, logger):
    logger.info("Preparing SNP INFO using NGT....")
    snp_info_by_quality = None
    try:
        mb_vcf.nonref_cells_by_numerical_genotype(task_result)
        mb_vcf.snp_info()
        mb_vcf.map_snp_info_to_amplicon(amp_bedfile)
        mb_vcf.add_ids_to_snp_info()
        snp_info_by_quality = mb_vcf.filter_snp_info_by_quality()
    except Exception as e:
        msg = "Failed with {} in preapring SNP INFO using NGT".format(
            str(e))
        logger.error(msg)

    if snp_info_by_quality is None:
        msg = "SNP_INFO shouldn't be None. It seems that \
            an wrong metadata.csv is given."
        raise ValueError(msg)

    logger.info("SNP INFO is generated successfully.")
    return snp_info_by_quality


def terminate_processes(proc_queue, logger, terminate=False):
    logger.info("Joining child processes......")
    for proc in proc_queue:
        logger.info("process {} ends".format(proc.pid))
        if terminate:
            proc.terminate()
        else:
            proc.join()
    logger.info("All child processes are joined sucessfully.")


def make_cli_parser():
    """Create the command line interface.

    :returnL an :py:class: 'argparse.ArgumentParser' instance

    """
    cli_parser = argparse.ArgumentParser(description=__doc__)
    cli_parser.add_argument(
        '--logfile', default='vcf_to_loom.log',
        help='Path to a log file'
    )
    cli_parser.add_argument(
        '--quality_threshold', type=int, default=1000,
        help='Quality threshold to seperate a loom file logically.'
    )
    cli_parser.add_argument(
        '--min_cell_mutated', type=int, default=5,
        help='mininum number of cells allowed to be mutated..'
    )
    cli_parser.add_argument(
        '--min_cell_percentage', type=float, default=0.001,
        help='Percentage of cells allowed to be mutated..'
    )
    cli_parser.add_argument(
        '--multi_processing', type=bool, default=False,
        help='A flag to run the module in parallel. (DEPRECATED)'
    )
    cli_parser.add_argument(
        'dir_yaml_file', type=str,
        help='Path to DIR.{genome}.yaml'
    )
    cli_parser.add_argument(
        'amplicon_bed_file', help='Path to an amplicon BED file'
    )
    cli_parser.add_argument(
        'metadata_csv_file', help='Path to an metadata CSV file'
    )
    cli_parser.add_argument(
        'vcf_file', help='Path to an VCF file'
    )
    cli_parser.add_argument(
        'output_filename', help='Path to an output file'
    )
    return cli_parser


def main(argv=None):
    cli_parser = make_cli_parser()
    args = cli_parser.parse_args(argv)
    configure_root_logger(filename=args.logfile)

    quality_threshold = args.quality_threshold
    # Disable multiprocessing: Causes sporadic failure where underlying data would be garbage collected and
    #  causes "Boolean index not matching dimension" with garbage length of array (random number)
    multi_processing = False
    dir_yaml = args.dir_yaml_file
    amp_bedfile = args.amplicon_bed_file
    metadata_csv = args.metadata_csv_file
    vcf_file = args.vcf_file
    output_filename = args.output_filename

    logger.info("VCF to Loom converstion starts")
    logger.info("Quality Threshold: {}".format(quality_threshold))
    logger.info("DIR.yaml: {}".format(dir_yaml))
    logger.info("Amplicon file: {}".format(amp_bedfile))
    logger.info("metadata csv: {}".format(metadata_csv))
    logger.info("VCF file: {}".format(vcf_file))
    logger.info("Output filename: {}".format(output_filename))

    start_time = time.time()

    global mb_vcf
    global data_dict
    global threadLock

    mb_vcf = None
    data_dict = {}
    with io.FileIO(vcf_file) as vcf_fh:
        mb_vcf = MBVcf(vcf_fh, quality_threshold,
                       REQUIRED_FIELDS_IN_VCF)

    num_of_samples = mb_vcf.num_of_samples
    logger.info("Number of samples: {}".format(num_of_samples))

    threads = []
    try:
        logger.info("Starting threads......")
        threadLock = threading.Lock()
        for task, fh_type in zip(['NGT', 'DP', 'GQ', 'AD', 'RO'],
                                 [amp_bedfile, None, None, None, None]):
            if multi_processing:
                tid = threading.Thread(target=thread_task,
                                       args=(task, fh_type, logger))
                tid.start()
                threads.append(tid)
                logger.info("Thread {} for {} start.".format(tid, task))
            else:
                thread_task(task, fh_type, logger)

        logger.info("Process(es) have been started.")

        logger.info("Waiting for threads to complete processing.....")
        for t in threads:
            t.join()
    except Exception as e:
        logger.info(str(e))
        raise Exception(str(e))

    snp_info = prepare_snp_info(data_dict['NGT'],
                                amp_bedfile, logger)
    data_dict['SNP_INFO'] = snp_info
    logger.info("All required data has been collected.")

    allsites_block = int(mb_vcf.allsites_block)
    logger.info("allsites_block starts at {}".format(allsites_block))

    logger.info("Reading genome version from YAML....")
    yaml_json = yaml.load(open(dir_yaml))

    genome_version = None
    try:
        genome_version = yaml_json['genome']['version']
    except Exception as e:
        msg = "Genome version does not exist in {}".format(dir_yaml)
        logger.info(msg)
        raise ValueError(msg)

    logger.info("Loading metadata CSV file.....")
    loom_meta_data = create_meta_data(amp_bedfile, metadata_csv,
                                      genome_version)
    loom_meta_data['allsites_block'] = allsites_block
    logger.info("metadata CSV has been loaded successfully.")

    col_attrs = {}
    barcodes = mb_vcf.barcodes
    if barcodes is not None:
        col_attrs["barcode"] = barcodes

    logger.info("Writing Loom file using data dict.....")
    LoomWriter.write(
        TapData(
            data_dict['DP'], data_dict['AD'], data_dict['GQ'],
            data_dict['NGT'], data_dict['RO'], data_dict['SNP_INFO']),
        output_filename,
        file_attrs=loom_meta_data,
        col_attrs=col_attrs
    )
    logger.info("Loom file is created successfully.")

    end_time = time.time()
    logger.info("---- Processing VCF file takes {}s seconds ----".format(
        end_time - start_time))
    logger.info("VCF to Loom conversion completes.")


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print("Exiting abnormally due to {}".format(str(e)))
