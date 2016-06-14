#!/usr/bin/env python
"""
Calculate distance to the boundary based on provided
coordinate file
"""

__date__ = "2014-12-19"
__author__ = "Rafal Gumienny"
__email__ = "r.gumienny@unibas.ch"
__license__ = "GPL"

# imports
import re
import sys
import time
import gzip
import subprocess
from sys import exit
from Bio import SeqIO
from os.path import abspath
from argparse import ArgumentParser, RawTextHelpFormatter

parser = ArgumentParser(description=__doc__, formatter_class=RawTextHelpFormatter)
parser.add_argument("-v",
                    "--verbose",
                    dest="verbose",
                    action="store_true",
                    default=False,
                    help="Be loud!")
parser.add_argument('--seq', default='seqs.fa', dest='seq',
                  action='store', help='fasta with mRNA sequences')
parser.add_argument('--out', '-o', default='output.tab', dest='out',
                  action='store', help='output table')
parser.add_argument('--coords', default='coords.tab', dest='coords',
                  action='store', help='file with target\
                  sites positions, miRNA and target gene ID')
parser.add_argument('--contextLen', type=int, dest='contextLen', default=50,
                     action='store', help='length of the \
        context sequence serounding binding site')


# redefine a functions for writing to stdout and stderr to save some writting
syserr = sys.stderr.write
sysout = sys.stdout.write


def main(options):
    """Main logic of the script"""
    coords = []
    if options.verbose == True:
            print "Reading sequences from coordinate file %s" % (options.coords)
    try:
        corfile = gzip.open(options.coords, 'rb')
    except IOError:
        raise IOError('Cannot read from coordinate file %s' % (options.coords))

    for c in corfile:
        try:
            c = c.rstrip().split()
            coords.append([c[0], int(c[2]), int(c[3]), c[1]])
        except ValueError, e:
            raise ValueError("Wrong coordinates: %s" % " ".join(c))


    # Read mRNA sequences into hash table: {id:sequence}
    if options.verbose == True:
            print "Reading sequences from mRNA file %s" % (options.seq)
    try:
        seq_obj = open(options.seq, 'Ur')
    except IOError:
        raise IOError('Cannot read from mRNA file %s' % (options.seq))

    mRNAseqs = {}
    for seq in SeqIO.parse(seq_obj, 'fasta'):
        mRNAseqs[str(seq.id)] = str(seq.seq)


    # Open output file and write first lines
    try:
        outfile = gzip.open(options.out, 'wb')
    except IOError:
        raise IOError("Connot open output file %s" % (options.out))
    # outfile.write('#siteID\tdistToBoundary\n')

    # Iterate through the binding coordinates to calculate their score
    if options.verbose == True:
            print "Calculating average G content... "
    for (mrnaid, lowerix, upperix, mirnas) in coords:
        if mrnaid in mRNAseqs:
            mrnasequ = mRNAseqs[mrnaid]
            score = calculate_distance_to_boundary(len(mrnasequ), lowerix)
            for mirna in mirnas.split(","):
                outtext = '%s,%s,%i,%i\t%s\n' % (mrnaid,
                                                 mirna,
                                                 lowerix,
                                                 upperix,
                                                 score)
                outfile.write(outtext)
        else:
            for mirna in mirnas.split(","):
                outtext = '%s,%s,%i,%i\t%s\n' % (mrnaid,
                                                 mirna,
                                                 lowerix,
                                                 upperix,
                                                 "NA")
                outfile.write(outtext)


    outfile.close()



def calculate_distance_to_boundary(mRNAlen, begpos):
        return str(begpos if begpos < (mRNAlen - begpos) else (mRNAlen - begpos))

if __name__ == '__main__':
    try:
        try:
            options = parser.parse_args()
        except Exception, e:
            parser.print_help()
            sys.exit()
        if options.verbose:
            start_time = time.time()
            start_date = time.strftime("%d-%m-%Y at %H:%M:%S")
            syserr("############## Started script on %s ##############\n" %
                   start_date)
        main(options)
        if options.verbose:
            syserr("### Successfully finished in %i seconds, on %s ###\n" %
                   (time.time() - start_time,
                    time.strftime("%d-%m-%Y at %H:%M:%S")))
    except KeyboardInterrupt:
        syserr("Interrupted by user after %i seconds!\n" %
               (time.time() - start_time))
        sys.exit(-1)
