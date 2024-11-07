#!/usr/bin/env python3
import pysam

"""
   This script as a lighter version of coverageMaker.py - it accepts the output of samtools collate
   and does not expect to see any secondary alignments. Will just detect multiple primary alignments and sort
"""

def parse_inputs():
    """Parse command-line arguments.

    Returns:
        parser(:obj:`argparse.ArgumentParser`)
        args(:obj:`ArgumentParser.parse_args`): A namespace containing the script's input arguments.
    """
    parser = argparse.ArgumentParser(description='Classify reads as host, graft, both, neither, or ambiguous.')
    parser.add_argument('-b', '--bam', help='Input BAM or CRAM file with BWA alignments', required=True)
    parser.add_argument('-o', '--output-dir', help='Output directory', required=True)
    parser.add_argument('-q', '--mapq', help="Minimum mapping quality", default="30")
    parser.add_argument('-i', '--increment', help="Increment to shift pid, useful for cluster jobs", default="0")
    args = parser.parse_args()
    return parser, args


def sort_by_pos(original_bam, in_bam):
    """

    A routine which sorts a BAM file. This is the final task before running coverage analysis

    Args:
         original_bam(str): path to the original input bam to retrieve its basename
         in_bam(str): input bam

    """
    vetted_name = os.path.basename(original_bam)
    vetted_name = vetted_name.rstrip(".bam").rstrip(".cram")
    sort_by_pos_bam = '{:s}/{:s}_cleaned.bam'.format(output_dir, vetted_name)
    sort("-o", sort_by_pos_bam, in_bam)
    return sort_by_pos_bam


def write_log(out_dir, pid, message):
    logName = "{:s}/{:d}_processing.log".format(out_dir, pid)
    with open(logName, "a") as logf:
        logf.write(message + "\n")
    logf.close()


def initialize_file_input(arguments):
    """Set file-relate inputs to variables with concise names.

    Args:
        arguments (:obj:`ArgumentParser.parse_args`): A namespace containing the script's input arguments.
    Returns:
        input_bam(:obj:`pysam.AlignmentFile`): A bam file object of all reads aligned.
        out_dir(str): A directory for writing files.
        min_q(int); Minimum mapping quality allowed.

    """
    input_bam = arguments.bam
    out_dir = arguments.output_dir
    min_q = arguments.mapq
    inc = arguments.increment
    return input_bam, out_dir, int(min_q), int(inc)


def _my_valid_reads(bam, pid):
    """
    We assume that pairs of alignments for the same read will be far from each other
    when bam is sorted by coordinate
    Args:
       bam(str): input bam file
       pid(int): process id
    """
    seen_reads = dict(primary=[])
    current_read = None
    skipped_reads = 0

    for read in bam.fetch(until_eof=True):
        if current_read is not None and read.query_name != current_read:
            if len(seen_reads['primary']) > 0:
                for r in seen_reads['primary']:
                    yield r
                seen_reads['primary'].clear()
        current_read = read.query_name
        if not read.is_secondary and len(seen_reads['primary']) < 2 and read.is_proper_pair and not read.is_unmapped:
            seen_reads['primary'].append(read)
        else:
            skipped_reads += 1
    write_log(output_dir, pid, f"Reads skipped due to Merging: {skipped_reads: 12}")


def merge_bam(sorted_alignments, out_dir, pid):
    """
    This function uses generator to make sure we have alignments (primary or secondary) for unique reads

    Args:
      sorted_alignments: merged sorted alignments with uniq read sequences, alignments primary or secondary
      out_dir(str): output directory
      pid(int): process id

    """
    s_bam = AlignmentFile(sorted_alignments, "rb")
    m_bam = AlignmentFile("{:s}/{:d}_merged.bam".format(out_dir, pid), "wb", template=s_bam)
    for read in _my_valid_reads(s_bam, pid):
        m_bam.write(read)

    s_bam.close()
    m_bam.close()


def remove_file(file_path: str):
    try:
        os.remove(file_path)
    except OSError:
        pass


if __name__ == '__main__':
    import argparse
    import os
    from pysam import sort
    from pysam import AlignmentFile

    parser, args = parse_inputs()

    # Load a bam file and split
    p_id = os.getpid()
    input_bam_path, output_dir, minimum_q, increment = initialize_file_input(args)

    # priority-merge primary alignments
    merge_bam(input_bam_path, output_dir, p_id + increment)
    bam_merged = "{:s}/{:d}_merged.bam".format(output_dir, p_id + increment)

    bam_sorted_pos = sort_by_pos(input_bam_path, bam_merged)
    remove_file(bam_merged)
    finalMetrics = "\n" + pysam.flagstat(bam_sorted_pos)
    write_log(output_dir, p_id + increment, finalMetrics)
