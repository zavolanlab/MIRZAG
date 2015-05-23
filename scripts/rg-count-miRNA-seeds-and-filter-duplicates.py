#!/usr/bin/env python

"""
Count miRNA seed in the provided sequences according to chosen definition
"""
__author__ = "Rafal Gumienny"
__email__ = "r.gumienny@unibas.ch"
__license__ = "GPL"

# imports
import sys
import time
import re
import gzip
import pandas as pd
from argparse import ArgumentParser
from Bio import SeqIO, Seq

parser = ArgumentParser(description=__doc__)
parser.add_argument("-v",
                    "--verbose",
                    dest = "verbose",
                    action="store_true",
                    default=False,
                    help="Be loud!")
parser.add_argument("--motifs",
                    help="miRNA/siRNA sequences to use when scanning")
parser.add_argument("--seqs",
                    help="Sequences for scaning eg. 3' UTRs")
parser.add_argument("--how",
                    choices=("ElMMo", "TargetScan", "6-mer"),
                    default="TargetScan",
                    help="What definition for seed to use")
parser.add_argument("--output",
                    dest="output",
                    default="coords.tab",
                    help="Name of the output file , defaults to coords.tab")
parser.add_argument("--context",
                    dest="context",
                    type=int,
                    default=50,
                    help="Context for sequence to print, defaults to 50")
parser.add_argument("--split-by",
                    dest="split_by",
                    required=True,
                    help="Split id by the string")
parser.add_argument("--index-after-split",
                    dest="index_after_split",
                    type=int,
                    required=True,
                    help="After split take this column as new id, 0 based")

try:
    options = parser.parse_args()
except Exception, e:
    parser.print_help()

# redefine a functions for writing to stdout and stderr to save some writting
syserr = sys.stderr.write
sysout = sys.stdout.write

def main():
    """Main logic of the script"""
    motifs = {str(rec.id):Seq.Seq(str(rec.seq).upper().replace('U','T')) for rec in SeqIO.parse(options.motifs, 'fasta')}
    seqs = {str(rec.id):str(rec.seq).upper().replace('U','T') for rec in SeqIO.parse(options.seqs, 'fasta')}

    data = []
    for mirid, mir in motifs.iteritems():
        if options.how == "ElMMo":
            # ElMMo definition of the seed
            regex = r"%s|%s" %(str(motifs[mirid][0:7].reverse_complement()),str(motifs[mirid][1:8].reverse_complement()))
            for seqid, seq in seqs.iteritems():
                for match in re.finditer(regex, seq):
                    dist_to_boundary = calculate_distance_to_boundary(len(seq), match.span()[0])
                    if dist_to_boundary < 50:
                        continue
                    lowerIndex, upperIndex = get_indices(options.context,
                            match.span()[0],
                            match.span()[1])
                    if upperIndex >= len(seq):
                        mRNA = seq[-1 * options.context:]
                    elif lowerIndex < 0:
                        mRNA = seq[:options.context]
                    else:
                        mRNA = seq[lowerIndex:upperIndex]
                    if len(mRNA) == options.context:
                        data.append((seqid,
                                     mirid,
                                     match.span()[0],
                                     match.span()[1],
                                     mRNA))
                    else:
                        continue
        if options.how == "TargetScan":
            # TargetScan definition of the seed
            regex = r"%s|%s" %(str(motifs[mirid][1:8].reverse_complement()), str(motifs[mirid][1:7].reverse_complement()) + 'A')
            for seqid, seq in seqs.iteritems():
                for match in re.finditer(regex, seq):
                    dist_to_boundary = calculate_distance_to_boundary(len(seq), match.span()[0])
                    if dist_to_boundary < 50:
                        continue
                    lowerIndex, upperIndex = get_indices(options.context,
                                                         match.span()[0],
                                                         match.span()[1])
                    if upperIndex >= len(seq):
                        mRNA = seq[-1 * options.context:]
                    elif lowerIndex < 0:
                        mRNA = seq[:options.context]
                    else:
                        mRNA = seq[lowerIndex:upperIndex]
                    if len(mRNA) == options.context:
                        data.append((seqid,
                                     mirid,
                                     match.span()[0],
                                     match.span()[1],
                                     mRNA))
                    else:
                        continue
        if options.how == "6-mer":
            # 6-mer
            regex = r"%s" %(str(motifs[mirid][1:7].reverse_complement()))
            for seqid, seq in seqs.iteritems():
                for match in re.finditer(regex, seq):
                    dist_to_boundary = calculate_distance_to_boundary(len(seq), match.span()[0])
                    if dist_to_boundary < 50:
                        continue
                    lowerIndex, upperIndex = get_indices(options.context,
                                                         match.span()[0],
                                                         match.span()[1])
                    if upperIndex >= len(seq):
                        mRNA = seq[-1 * options.context:]
                    elif lowerIndex < 0:
                        mRNA = seq[:options.context]
                    else:
                        mRNA = seq[lowerIndex:upperIndex]
                    if len(mRNA) == options.context:
                        data.append((seqid,
                                     mirid,
                                     match.span()[0],
                                     match.span()[1],
                                     mRNA))
                    else:
                        continue

    names = ['id', 'mirna', 'beg', 'end', 'seq']
    data = pd.DataFrame(data, columns=names)
    data['newid'] = [i.split(options.split_by)[options.index_after_split] for i in data.id]
    ndf = data.drop_duplicates(cols=['newid', 'mirna', 'seq'])
    with gzip.open(options.output, 'wb') as o:
        ndf[names].to_csv(o, header=False, index=False, sep='\t')
    if options.verbose:
        syserr("Filtering:\n")
        syserr(" - number of coordinates before: %i\n" % len(data.index))
        syserr(" - number of coordinates after : %i\n" % len(ndf.index))


def get_indices(desired_length, beg, end):
    """Calculate upper and lower index based on what size we want and what are the coordinates"""
    lower_index = beg - ((desired_length - (end - beg))//2 + (desired_length - (end - beg))%2)
    upper_index = end + ((desired_length - (end - beg))//2)
    return lower_index, upper_index

def calculate_distance_to_boundary(mRNAlen, begpos):
    "Calculate distance to the closest boundary"
    return str(begpos if begpos < (mRNAlen - begpos) else (mRNAlen - begpos))

if __name__ == '__main__':
    start_time = time.time()
    try:
        if options.verbose:
            start_date = time.strftime("%d-%m-%Y at %H:%M:%S")
            syserr("############## Started script on %s ##############\n" % start_date)
        main()
        if options.verbose:
            syserr("### Successfully finished in %i seconds, on %s ###\n" % (time.time() - start_time, time.strftime("%d-%m-%Y at %H:%M:%S")))
    except KeyboardInterrupt:
        syserr("Interrupted by user after %i seconds!\n" % (time.time() - start_time))
        sys.exit(-1)
