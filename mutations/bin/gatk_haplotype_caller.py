#!/usr/bin/env python3

import os
import re
import subprocess
import multiprocessing
from functools import partial
from pathlib import Path

import click

from missionbio.dna.config import DEFAULTS, Config
from missionbio.dna.preflight import validate_panel

HOST_ID_LOC = 2

BARCODE_RE = re.compile(r"@RG\s+ID:([^\s]+)")


def convert_to_bam_with_index(sam_file, bam_file, delete=True):
    """Convert sam file to bam file

    Parameters
    ----------
    sam_file: Path
        path to the file to convert
    bam_file: Path
        path where the bam file will be created
    delete: bool
        when True, sam file is removed after conversion
    """
    if not sam_file.exists():
        raise ValueError(f"Sam file {sam_file} does not exist")

    subprocess.run(["samtools", "view", "-h", "-b", sam_file, "-o", bam_file], check=True)
    subprocess.run(["samtools", "index", bam_file])
    if delete:
        sam_file.unlink()


def get_barcode(sam_file):
    """Get barcode from sam_file header

    Raises ValueError if sam_file contains more than one barcode

    Parameters
    ----------
    sam_file: Path
        path to the sam file

    Returns
    -------
    str
        barcode associated with this file
    """
    result = subprocess.run(
        ["samtools", "view", "-H", sam_file], check=True, stdout=subprocess.PIPE, encoding="utf-8"
    )
    barcodes = BARCODE_RE.findall(result.stdout)
    if len(barcodes) != 1:
        how_many = "more than one" if len(barcodes) > 1 else "no"
        raise ValueError(f"Sam file {sam_file} contains {how_many} variant\n({barcodes[:5]})")

    return barcodes[0]


def read_gatk_params(config):
    """Read gatk params from run options file

    Parameters
    ----------
    config: Dict
        path to the config.yaml

    Returns
    -------
    List[str]
        list of gatk command line options
    """
    gatk_params = []

    if "gatk" in config:
        for group in ["common", "HaplotypeCaller"]:
            for key, value in config["gatk"].get(group, {}).items():
                if value:
                    gatk_params.append(f"--{key}")
                    gatk_params.append(f"{value}")

    return gatk_params


@click.command()
@click.option(
    "--input-dir",
    type=click.Path(exists=True, file_okay=False),
    help="Path to the dir with .sam files",
)
@click.option(
    "--input-file", type=click.Path(exists=True, dir_okay=False), help="Path to a single .sam file"
)
@click.option(
    "--output-dir",
    required=True,
    type=click.Path(file_okay=False),
    help="Path where the GVCF files will be created",
)
@click.option(
    "--cores",
    "cores",
    required=True,
    help="The number of cores used",
)
@click.option(
    "--config",
    "config_file",
    required=True,
    type=click.Path(exists=True, dir_okay=False),
    help="Path to the DIR.yaml file",
)
@click.option("--prefix", help="When set, only files starting with prefix will be processed")
def main(input_dir, input_file, output_dir, cores, config_file, prefix):
    if int(os.environ.get("DISABLE_GENOTYPE_CALLING", 0)):
        print("GATK GENOTYPE CALLING is disabled.")
        return

    if input_dir is None:
        if input_file is None:
            raise click.UsageError("Please specify --input-dir or --input-file")
        else:
            inputs = [Path(input_file).absolute()]
    else:
        if input_file is None:
            print(f"Reading .sam files with prefix {prefix} from {input_dir}")
            inputs = [f.absolute() for f in Path(input_dir).glob(f"{prefix}*.sam")]
        else:
            raise click.UsageError(
                "Only one of --input-dir, --input-file can be used at the same time"
            )

    output_dir = Path(output_dir)
    config_file = Path(config_file)
    config = Config.read(config_file, DEFAULTS)
    validate_panel(config)
    reference = config["genome"]["path"]
    intervals = config["panel"]["bed"]
    gatk_params = read_gatk_params(config)

    output_dir.mkdir(parents=True, exist_ok=True)

    pool = multiprocessing.Pool(int(cores))
    pool.map(partial(run_gatk, reference=reference, intervals=intervals, output_dir=output_dir, gatk_params=gatk_params), inputs)


def run_gatk(input_sam, reference, intervals, output_dir, gatk_params):
    barcode = get_barcode(input_sam)
    bam_file = input_sam.parent / f"{barcode}.bam"
    gvcf_file = output_dir / f"{barcode}.g.vcf.gz"
    convert_to_bam_with_index(input_sam, bam_file)

    subprocess.run(
        [
            "gatk",
            "--java-options",
            "-Xmx2g",
            "HaplotypeCaller",
            "-R",
            reference,
            "-I",
            bam_file,
            "--emit-ref-confidence",
            "BP_RESOLUTION",
            "--native-pair-hmm-threads",
            "2",  # multiple instances will run in parallel
            "-L",
            intervals,
            "--verbosity",
            "ERROR",
            "-O",
            gvcf_file,
        ]
        + gatk_params,
        check=True,
    )
    print(f"GATK HC for {barcode} completed successfully")

if __name__ == "__main__":
    main()
