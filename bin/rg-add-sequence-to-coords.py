#!/usr/bin/env python
"""

"""

__date_ = "2014-08-05"
__author__ = "Rafal Gumienny"
__email__ = "r.gumienny@unibas.ch"
__license__ = "GPL"

# imports
import sys
import csv
import time
from Bio import SeqIO
from argparse import ArgumentParser

parser = ArgumentParser(description=__doc__)
parser.add_argument("-v",
                    "--verbose",
                    dest="verbose",
                    action="store_true",
                    default=False,
                    help="Be loud!")
parser.add_argument("--coords",
                    dest="coords",
                    required=True,
                    help="File with coordinates")
parser.add_argument("--seqs",
                    dest="seqs",
                    required=True,
                    help="UTR sequences used to generate coordinate file")
parser.add_argument("--output",
                    dest="output",
                    default="coords.tab",
                    help="Name of the output file , defaults to coords.tab")
parser.add_argument("--context",
                    dest="context",
                    type=int,
                    default=50,
                    help="Context for sequence to print, defaults to 50")

try:
    options = parser.parse_args()
except Exception, e:
    parser.print_help()

# redefine a functions for writing to stdout and stderr to save some writting
syserr = sys.stderr.write
sysout = sys.stdout.write


def main():
    """Main logic of the script"""
    mRNAseqs = read_fasta_to_dict(options.seqs)
    if options.verbose:
        syserr("Adjusting data\n")
    all_count = 0
    valid_count = 0
    with open(options.coords) as f, open(options.output, 'w') as out:
        for seqid, mirid, beg, end in csv.reader(f, delimiter='\t'):
            beg = int(beg)
            end = int(end)
            all_count += 1
            seq = mRNAseqs[seqid]
            dist_to_boundary = calculate_distance_to_boundary(len(seq), beg)
            if dist_to_boundary < 50:
                continue
            lowerIndex, upperIndex = get_indices(options.context,
                                                 beg,
                                                 end)
            if upperIndex >= len(seq):
                mRNA = seq[-1 * options.context:]
            elif lowerIndex < 0:
                mRNA = seq[:options.context]
            else:
                mRNA = seq[lowerIndex:upperIndex]
            if len(mRNA) == options.context:
                out.write("%s\t%s\t%i\t%i\t%s\n" % (seqid,
                                                    mirid,
                                                    beg,
                                                    end,
                                                    mRNA))
            valid_count += 1
    if options.verbose:
        syserr("Wrote %i coordinates out of %i\n" % (valid_count, all_count))


def get_indices(desired_length, beg, end):
    """Calculate upper and lower index based on what size we want and what are the coordinates"""
    lower_index = beg - ((desired_length - (end - beg))//2 + (desired_length - (end - beg))%2)
    upper_index = end + ((desired_length - (end - beg))//2)
    return lower_index, upper_index


def calculate_distance_to_boundary(mRNAlen, begpos):
    "Calculate distance to the closest boundary"
    return begpos if begpos < (mRNAlen - begpos) else (mRNAlen - begpos)


def read_fasta_to_dict(path_to_file):
    """Read fasta file into dictionary

    Args:
        path_to_file (@todo): @todo

    Returns: @todo

    """
    if options.verbose:
        syserr("Reading sequences from %s \n" % (path_to_file))
    try:
        seq_obj = open(path_to_file, 'Ur')
        seqs = {}
        for seq in SeqIO.parse(seq_obj, 'fasta'):
            seqs[str(seq.id)] = str(seq.seq)
    except IOError:
        raise IOError('Cannot read from %s' % (path_to_file))

    return seqs

if __name__ == '__main__':
    try:
        if options.verbose:
            start_time = time.time()
            start_date = time.strftime("%d-%m-%Y at %H:%M:%S")
            syserr("############## Started script on %s ##############\n" % start_date)
        main()
        if options.verbose:
            syserr("### Successfully finished in %i seconds, on %s ###\n" % (time.time() - start_time, time.strftime("%d-%m-%Y at %H:%M:%S")))
    except KeyboardInterrupt:
        syserr("Interrupted by user after %i seconds!\n" % (time.time() - start_time))
        sys.exit(-1)
