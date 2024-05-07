#!/usr/bin/env python3
import multiprocessing
import sys
from argparse import ArgumentParser, Namespace
from pathlib import Path
from subprocess import PIPE, CalledProcessError, Popen


def main(args):
    print("splitBAMByBarcode is starting......")

    # Create a queue and worker for flushing output files
    write_queue = multiprocessing.JoinableQueue(args.max_pending_writes)
    worker = Worker(write_queue)
    worker.start()

    try:
        writers = {}
        header = read_header(args.input)
        for barcode, line in read_raw_barcodes(args.input):
            if barcode not in writers:
                writers[barcode] = SamWriter(
                    barcode,
                    args.output_dir,
                    args.num_of_gatk,
                    write_queue=write_queue,
                    batch_size=args.batch_size,
                )
                for header_line in header:
                    if header_line.startswith("@RG"):
                        if barcode not in header_line:
                            continue
                    writers[barcode].write(header_line)
            writers[barcode].write(line)
        print("Flushing output files")
        for w in writers.values():
            w.flush(final=True)
    finally:
        # Send poison pill to worker and wait for it to stop
        write_queue.put(None)
        write_queue.join()

    print("splitBAMByBarcode completes successfully.")


def read_header(filename):
    """Read bam header

    Parameters
    ----------
    filename: str
        file to read the header from

    Returns
    -------
    List[str]
        list of header lines
    """
    cmd = " ".join(["samtools", "view", "-H", filename])
    header = []
    with Popen(cmd, stdout=PIPE, shell=True, encoding="utf-8") as proc:
        for index, line in enumerate(proc.stdout):
            header.append(line.strip())
        if proc.wait() != 0:
            raise CalledProcessError(proc.returncode, cmd)
    return header


def read_raw_barcodes(filename):
    """Yield all raw barcodes found in file

    Equivalent to
    ```
    sam_fh = pysam.AlignmentFile(samfile, mode='rb')
    for index, alignment_segment in enumerate(sam_fh.fetch(until_eof=True)):
        yield alignment_segment.get_tag('RG'), index
    ```
    but about 50% faster

    Starts a separate process that extracts barcodes. The
    barcodes are yielded as soon as the process prints them
    to the stdout.

    Parameters
    ----------
    filename: str
        path to the bam file

    Returns
    -------
    Generator[Tuple[str, str]]
    """
    find_RG = """\
    '$_ =~ /RG:Z:([^\\s]+)/;
    print "$1 $_"'
    """
    # The following command runs samtools to read alignments
    # and pipes them through a perl script which extracts unique barcodes
    cmd = " ".join(["samtools", "view", filename, "|", "perl", "-ne", find_RG])
    with Popen(cmd, stdout=PIPE, shell=True, encoding="utf-8") as proc:
        for index, line in enumerate(proc.stdout):
            if index % 100000 == 0:
                print(f"\rProcessed {index/1000000}M reads", end="", file=sys.stderr)
            barcode, _, line = line.strip().partition(" ")
            if not barcode:
                raise ValueError(f"Could not find barcode in read {line}")
            yield barcode, line
        print(f"\rProcessed {index/1000000}M reads", file=sys.stderr)
        if proc.wait() != 0:
            raise CalledProcessError(proc.returncode, cmd)


class SamWriter:
    """
    Batch writer for sam files

    Stores up to batch_size lines in memory and adds them to write_queue
    when bach_size is reached.
    """

    FILENAME_TEMPLATE = "{output_dir}/gatk_hc_{host}_{barcode}.sam"
    host = -1

    def __init__(
        self, barcode, output_dir=".", number_of_hosts=1, write_queue=None, batch_size=1000
    ):
        """
        Parameters
        ----------
        barcode: str
            barcode to create the file for
        output_dir: str
            directory where output files are created
        number_of_hosts: int
            number of gatk hosts
        write_queue: multiprocessing.JoinableQueue
            queue for scheduling writes
        batch_size: int
            number of indices to store in memory before caching to disk
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        self.filename = self.FILENAME_TEMPLATE.format(
            output_dir=output_dir, host=self.get_host_id(number_of_hosts), barcode=barcode
        )
        self.batch_size = batch_size
        self.first_batch = True
        self.lines = []
        with open(self.filename, "w"):
            pass
        self.write_queue = write_queue

    def write(self, index):
        self.lines.append(index)
        if len(self.lines) > self.batch_size:
            self.flush()

    def flush(self, final=False):
        self.write_queue.put((self.filename, self.lines, self.first_batch, final))
        self.first_batch = False
        self.lines = []

    @classmethod
    def get_host_id(cls, number_of_hosts):
        cls.host = (cls.host + 1) % number_of_hosts
        return cls.host + 1


class Worker(multiprocessing.Process):
    """
    Worker for flushing sam files
    """

    def __init__(self, write_queue):
        """
        Parameters
        ----------
        write_queue: multiprocessing.JoinableQueue
            queue to read write requests from
        """
        super().__init__()
        self.write_queue = write_queue

    def run(self):
        """
        Keep reading write requests from queue until poison pill (None)
        is encountered

        Write requests have the following structure:
        (filename, lines, first_batch, final)

        filename: str
            name of the file to write the indices to
        lines: List[str]
            lines to write to file
        first_batch: bool
            True if this is the first batch of lines
        final: bool
            True if this is the final write to the file
        """

        while True:
            next_task = self.write_queue.get()
            # check for poison pill
            if next_task is None:
                self.write_queue.task_done()
                break

            # decode write request, write changes to file
            filename, lines, first_batch, final = next_task
            with open(filename, "a") as f:
                for line in lines:
                    f.write(line)
                    f.write("\n")
            self.write_queue.task_done()


class Parameters(Namespace):
    input = None
    output_dir = None
    num_of_gatk = None
    batch_size = 1000
    max_pending_writes = 10000

    @staticmethod
    def from_args():
        parser = ArgumentParser(description="Split BAM By Barcode")
        parser.add_argument("--input", "-i", required=True, help="path to an input BAM file")
        parser.add_argument(
            "--output-dir", "-o", required=True, help="output directory for read group"
        )
        parser.add_argument(
            "--num-of-gatk", "-n", type=int, required=True, help="number of GATK instances"
        )
        parser.add_argument(
            "--batch-size",
            "-b",
            type=int,
            default=1000,
            help="Maximum number of reads to keep in memory before flushing to file",
        )
        parser.add_argument(
            "--max-pending-writes",
            type=int,
            default=Parameters.max_pending_writes,
            help="Number of unflushed writes to allow. Set to 0 for no limit.",
        )
        args = parser.parse_args(namespace=Parameters())
        return args


if __name__ == "__main__":
    args = Parameters.from_args()
    main(args)
