#!/usr/bin/env python
import optparse
from Bio import SeqIO
from sys import exit
import sys
import gzip
from os.path import abspath
import subprocess
import re

# parse arguments and store them in the variables
parser = optparse.OptionParser()
parser.add_option('-v',
                     '--verbose',
                     dest='verbose',
                     action="store_true",
                     default=False)
parser.add_option('--seq',
                     default='seqs.fa',
                     dest='seq',
                     action='store',
                     help='fasta with mRNA sequences')
parser.add_option('--out',
                     '-o',
                     default='output.tab',
                     dest='out',
                     action='store',
                     help='output table')
parser.add_option('--coords',
                     default='coords.tab',
                     dest='coords',
                     action='store',
                     help='file with target\
                     sites positions, miRNA and target gene ID')
parser.add_option('--contextLen',
                     type=int,
                     dest='contextLen',
                     default=50,
                     action='store',
                     help='length of the \
                     context sequence serounding binding site')



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


def main(arguments, verbose):
    # Read coords file into the table: [geneID, miR name, begining position,
    # end position]
    coords = []
    if verbose == True:
        print "Reading sequences from coordinate file %s" % (arguments.coords)
    try:
        corfile = gzip.open(arguments.coords, 'r')
    except IOError:
        raise IOError( 'Cannot read from coordinate file %s' % (arguments.coords))

    for c in corfile:
        try:
            c = c.rstrip().split()
            coords.append([c[0], int(c[1]), int(c[2]), c[3]])
        except ValueError, e:
            raise ValueError("Wrong coordinates: %s" % " ".join(c))


    # Read mRNA sequences into hash table: {id:sequence}
    if verbose == True:
        print "Reading sequences from mRNA file %s" % (arguments.seq)
    try:
        seq_obj = open(arguments.seq, 'Ur')
    except IOError:
        raise IOError( 'Cannot read from mRNA file %s' % (arguments.seq))

    mRNAseqs = {}
    for seq in SeqIO.parse(seq_obj, 'fasta'):
        mRNAseqs[str(seq.id)] = str(seq.seq)


    # Open output file and write first lines
    try:
        outfile = gzip.open(arguments.out, 'wb')
    except IOError:
        raise IOError( "Connot open output file %s" % (arguments.out))
    # outfile.write('#siteID\tflanksG\tflanksA\tflanksC\tflanksU\n')

    # Iterate through the binding coordinates to calculate their score
    if verbose == True:
        print "Calculating average G content... "
    clen = float(len(coords))
    for (mrnaid, lowerix, upperix, mirnas) in coords:
        # here we assume that the coordinates are given in 0-based for start and
        # 1-based for end
        if mrnaid in mRNAseqs:
            mrnasequ = mRNAseqs[mrnaid]
            coorlen = upperix - lowerix
            lower_seq = mrnasequ[lowerix - arguments.contextLen: lowerix]
            upper_seq = mrnasequ[upperix: upperix + arguments.contextLen]
            seq_len = (len(lower_seq) + len(upper_seq) + coorlen)
            # print len(lower_seq), len(upper_seq), upperix - lowerix, seq_len
            if lowerix - arguments.contextLen >= 0 and upperix <= len(mrnasequ) \
                and seq_len == (2 * arguments.contextLen + coorlen):
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


if __name__ == '__main__':
    arguments, args = parser.parse_args()
    verbose = arguments.verbose
    main(arguments, verbose)
