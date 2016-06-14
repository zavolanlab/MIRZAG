#!/usr/bin/env python
"""
Calculate flanks composition using provided coordinates
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
import optparse
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
parser.add_argument('--seq',
                  default='seqs.fa',
                  dest='seq',
                  action='store',
                  help='fasta with mRNA sequences')
parser.add_argument('--out',
                  default='output.tab',
                  dest='out',
                  action='store',
                  help='output table')
parser.add_argument('--coords',
                  default='coords.tab',
                  dest='coords',
                  action='store',
                  help='file with target\
                  sites positions, miRNA and target gene ID')
parser.add_argument('--contextLen',
                     type=int,
                     dest='contextLen',
                     default=50,
                     action='store',
                     help='length of the \
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
        corfile = gzip.open(options.coords, 'r')
    except IOError:
        raise IOError( 'Cannot read from coordinate file %s' % (options.coords))

    for c in corfile:
        try:
            c = c.rstrip().split()
            coords.append([c[0], int(c[2]), int(c[3]), c[1]])
        except ValueError:
            raise ValueError("Wrong coordinates: %s" % " ".join(c))


    # Read mRNA sequences into hash table: {id:sequence}
    if options.verbose == True:
        print "Reading sequences from mRNA file %s" % (options.seq)
    try:
        seq_obj = open(options.seq, 'Ur')
    except IOError:
        raise IOError( 'Cannot read from mRNA file %s' % (options.seq))

    mRNAseqs = {}
    for seq in SeqIO.parse(seq_obj, 'fasta'):
        mRNAseqs[str(seq.id)] = str(seq.seq)


    # Open output file and write first lines
    try:
        outfile = gzip.open(options.out, 'wb')
    except IOError:
        raise IOError( "Connot open output file %s" % (options.out))
    # outfile.write('#siteID\tflanksG\tflanksA\tflanksC\tflanksU\n')

    # Iterate through the binding coordinates to calculate their score
    if options.verbose == True:
        print "Calculating average G content... "
    clen = float(len(coords))
    for (mrnaid, lowerix, upperix, mirnas) in coords:
        # here we assume that the coordinates are given in 0-based for start and
        # 1-based for end
        if mrnaid in mRNAseqs:
            mrnasequ = mRNAseqs[mrnaid]
            coorlen = upperix - lowerix
            lower_seq = mrnasequ[lowerix - options.contextLen: lowerix]
            upper_seq = mrnasequ[upperix: upperix + options.contextLen]
            seq_len = (len(lower_seq) + len(upper_seq) + coorlen)
            # print len(lower_seq), len(upper_seq), upperix - lowerix, seq_len
            if lowerix - options.contextLen >= 0 and upperix <= len(mrnasequ) \
                and seq_len == (2 * options.contextLen + coorlen):
                gcont, acont, ccont, tcont = calculate_flanks_composition(lower_seq,
                                                                          upper_seq)
                for mirna in mirnas.split(","):
                    outtext = '%s,%s,%i,%i\t%s\t%s\t%s\t%s\n' % (mrnaid,
                                                                 mirna,
                                                                 lowerix,
                                                                 upperix,
                                                                 gcont,
                                                                 acont,
                                                                 ccont,
                                                                 tcont)
                    outfile.write(outtext)
            else:
                sys.stderr.write("Flanks out of boarders\n")
                for mirna in mirnas.split(","):
                    outtext = '%s,%s,%i,%i\t%s\t%s\t%s\t%s\n' % (mrnaid,
                                                                 mirna,
                                                                 lowerix,
                                                                 upperix,
                                                                 'NA',
                                                                 'NA',
                                                                 'NA',
                                                                 'NA')
                    outfile.write(outtext)
        else:
            for mirna in mirnas.split(","):
                outtext = '%s,%s,%i,%i\t%s\t%s\t%s\t%s\n' % (mrnaid,
                                                             mirna,
                                                             lowerix,
                                                             upperix,
                                                             'NA',
                                                             'NA',
                                                             'NA',
                                                             'NA')
                outfile.write(outtext)


    outfile.close()


def calculate_flanks_composition(lower_seq, upper_seq):
    extr_seq = lower_seq + upper_seq
    extr_seq = extr_seq.upper().replace('U', 'T')
    G_count = float(extr_seq.count('G'))
    A_count = float(extr_seq.count('A'))
    C_count = float(extr_seq.count('C'))
    T_count = float(extr_seq.count('T'))
    g_perc = str(G_count / float(len(extr_seq)))
    a_perc = str(A_count / float(len(extr_seq)))
    c_perc = str(C_count / float(len(extr_seq)))
    t_perc = str(T_count / float(len(extr_seq)))
    return g_perc, a_perc, c_perc, t_perc

if __name__ == '__main__':
    try:
        try:
            options = parser.parse_args()
        except Exception, e:
            parser.print_help()
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
