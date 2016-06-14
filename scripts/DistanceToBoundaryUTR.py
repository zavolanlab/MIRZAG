#!/usr/bin/env python
import optparse
import gzip
from Bio import SeqIO
from sys import exit
from os.path import abspath
import subprocess
import re

# parse arguments and store them in the variables
parser = optparse.OptionParser()
parser.add_option('-v', '--verbose', dest='verbose', action="store_true",
                     default=False)
parser.add_option('--seq', default='seqs.fa', dest='seq',
                     action='store', help='fasta with mRNA sequences')
parser.add_option('--out', '-o', default='output.tab', dest='out',
                     action='store', help='output table')
parser.add_option('--coords', default='coords.tab', dest='coords',
                     action='store', help='file with target\
        sites positions, miRNA and target gene ID')
parser.add_option('--contextLen', type=int, dest='contextLen', default=50,
                     action='store', help='length of the \
        context sequence serounding binding site')



def calculate_distance_to_boundary(mRNAlen, begpos):
        return str(begpos if begpos < (mRNAlen - begpos) else (mRNAlen - begpos))


def main(arguments, verbose):
    # Read coords file into the table: [geneID, miR name, begining position,
    # end position]
    coords = []
    if verbose == True:
            print "Reading sequences from coordinate file %s" % (arguments.coords)
    try:
        corfile = gzip.open(arguments.coords, 'rb')
    except IOError:
        raise IOError('Cannot read from coordinate file %s' % (arguments.coords))

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
        raise IOError('Cannot read from mRNA file %s' % (arguments.seq))

    mRNAseqs = {}
    for seq in SeqIO.parse(seq_obj, 'fasta'):
        mRNAseqs[str(seq.id)] = str(seq.seq)


    # Open output file and write first lines
    try:
        outfile = gzip.open(arguments.out, 'wb')
    except IOError:
        raise IOError("Connot open output file %s" % (arguments.out))
    # outfile.write('#siteID\tdistToBoundary\n')

    # Iterate through the binding coordinates to calculate their score
    if verbose == True:
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
if __name__ == '__main__':
    arguments, args = parser.parse_args()
    verbose = arguments.verbose
    main(arguments, verbose)
