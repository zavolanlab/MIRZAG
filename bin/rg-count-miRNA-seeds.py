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
parser.add_argument("--coords",
                    action="store_true",
                    default=False,
                    help="Print in simple format")
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


def get_indices(desired_length, beg, end):
    """Calculate upper and lower index based on what size we want and what are the coordinates"""
    lower_index = beg - ((desired_length - (end - beg))//2 + (desired_length - (end - beg))%2)
    upper_index = end + ((desired_length - (end - beg))//2)
    return lower_index, upper_index

def main():
    """Main logic of the script"""
    syserr("Reading motifs\n")
    motifs = {str(rec.id):Seq.Seq(str(rec.seq).upper().replace('U','T')) for rec in SeqIO.parse(options.motifs,'fasta')}
    syserr("Reading sequences\n")
    seqs = {str(rec.id):str(rec.seq).upper().replace('U','T') for rec in SeqIO.parse(options.seqs,'fasta')}
    seqslen = len(seqs.keys())

    seqcout = 1
    for seqid, seq in seqs.iteritems():
        syserr("Scanning %i sequence out of %i\r" % (seqcout, seqslen))
        seqcout += 1
        for mirid, mir in motifs.iteritems():
            if options.how == "ElMMo":
                # ElMMo definition of the seed
                regex = r"%s|%s" %(str(motifs[mirid][0:7].reverse_complement()),str(motifs[mirid][1:8].reverse_complement()))
            if options.how == "TargetScan":
                # TargetScan definition of the seed
                regex = r"%s|%s" %(str(motifs[mirid][1:8].reverse_complement()),str(motifs[mirid][1:7].reverse_complement()) + 'A')
            if options.how == "6-mer":
                # 6-mer
                regex = r"%s" %(str(motifs[mirid][1:7].reverse_complement()))
            for match in re.finditer(regex, seq):
                lowerIndex, upperIndex = get_indices(options.context,
                        match.span()[0],
                        match.span()[1])
                if upperIndex >= len(seq):
                    mRNA = seq[-1 * options.context:]
                elif lowerIndex < 0:
                    mRNA = seq[:options.context]
                else:
                    mRNA = seq[lowerIndex:upperIndex]

                if options.coords:
                    sysout("%s\t%s\t%i\t%i\t%s\n" % (seqid, mirid,
                        match.span()[0], match.span()[1], mRNA))
                else:
                    sysout("%s\t%s\t%i\t%i\t%s\t%s\n" % (seqid, mirid, match.span()[0], match.span()[1], regex, match.group()))
    syserr("\n")

if __name__ == '__main__':
    try:
        start_time = time.time()
        start_date = time.strftime("%d-%m-%Y at %H:%M:%S")
        syserr("############## Started script on %s ##############\n" % start_date)
        main()
        syserr("### Successfully finished in %i seconds, on %s ###\n" % (time.time() - start_time, time.strftime("%d-%m-%Y at %H:%M:%S")))
    except KeyboardInterrupt:
        syserr("Interrupted by user after %i seconds!\n" % (time.time() - start_time))
        sys.exit(-1)
