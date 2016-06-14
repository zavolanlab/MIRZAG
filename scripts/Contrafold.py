#!/usr/bin/env python
import os
import gzip
import shutil
import optparse
import subprocess
from Bio import SeqIO
from itertools import izip

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
parser.add_option('--contextLen_L', type=int, dest='contextLen_L',
                     default=0,
                     action='store', help='length of the \
                     context sequence downstream binding site to be unwinded')
parser.add_option('--contextLen_U', type=int, dest='contextLen_U',
                     default=0,
                     action='store', help='length of the \
                     context sequence upstream binding site to be unwinded')
parser.add_option('--context', type=int, dest='context',
                     default=50,
                     action='store', help='length of the \
                     context of the seed to be checked')
parser.add_option('--contrabin',
                     default="contrafold",
                     help="Path to CONTRAfold binary")


def main(arguments, verbose):
    # Read coords file into the table: [geneID, miR name, begining position,
    # end position]
    coords = []
    if verbose is True:
        print "Reading sequences from coordinate file %s" % (arguments.coords)
    try:
        corfile = gzip.open(arguments.coords)
    except IOError:
        raise IOError('Cannot read from coordinate file %s' % (arguments.coords))

    for c in corfile:
        try:
            c = c.rstrip().split()
            coords.append([c[0], int(c[1]), int(c[2]), c[3]])
        except ValueError, e:
            raise ValueError("Wrong coordinates: %s" % " ".join(c))


    # Read mRNA sequences into hash table: {id:sequence}
    if verbose is True:
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
        raise IOError( "Connot open output file %s" % (arguments.out))
    # outfile.write('#siteID\tContraScoreTargetSite\n')

    # Iterate through the binding coordinates to calculate their score
    if verbose is True:
        print "Calculating Contrafold accessibility..."
    count = 0
    corlen = len(coords)

    # make a directory to store the files
    dirpath = os.path.splitext(os.path.abspath(arguments.coords))[0]
    os.mkdir(dirpath)
    bad_coords = []
    mirnas_dict = {}
    for (mrnaid, lowerix, upperix, mirnas) in coords:
        id_to_write = ",".join(str(c) for c in [mrnaid, lowerix, upperix])
        inputname = os.path.join(dirpath, id_to_write + ".bpseq")
        if id_to_write in mirnas_dict:
            raise Exception("%s already in database" % id_to_write)
        mirnas_dict[id_to_write] = mirnas.split(",")
        # we assume that coordinates are like in bed file: start 0-based and end 1-based
        if mrnaid in mRNAseqs:
            mrnasequ = mRNAseqs[mrnaid]
            coorLen = upperix - lowerix
            lowerIndex = arguments.context - arguments.contextLen_L
            upperIndex = arguments.context + arguments.contextLen_U + coorLen
            mrna_fragment = mrnasequ[lowerix - arguments.context:upperix + arguments.context]
            mrnarange = range(1, len(mrna_fragment) + 1)
            if len(mrnasequ) >= upperIndex and  lowerIndex >= 0 \
                    and len(mrna_fragment) == 2 * arguments.context + coorLen:
                with open(inputname, 'w') as mrnainput:
                    bpseq = []
                    try:
                        for i in mrnarange:
                            if i >= lowerIndex + 1 and i <= upperIndex:
                                bpseq.append('%i\t%s\t%i' % (i, mrna_fragment[i - 1], 0))
                            else:
                                bpseq.append('%i\t%s\t%i' % (i, mrna_fragment[i - 1], -1))
                        mrnainput.write("\n".join(bpseq))
                    except IOError:
                        raise IOError('Cannot write into the MIRZA input files')
            else:
                bad_coords.append(id_to_write)
                continue
    contra_command_unwind = "%s predict %s --constraints --partition" % \
                            (arguments.contrabin,
                             os.path.join(dirpath, "*.bpseq"))
    contrarun_unwind = subprocess.Popen(contra_command_unwind,
                                        shell=True,
                                        stdout=subprocess.PIPE,
                                        stderr=subprocess.PIPE)
    stdout_unwind, stderr_unwind = contrarun_unwind.communicate()

    contra_command_wind = "%s predict %s --partition" % \
                          (arguments.contrabin,
                           os.path.join(dirpath, "*.bpseq"))
    contrarun_wind = subprocess.Popen(contra_command_wind,
                                      shell=True,
                                      stdout=subprocess.PIPE,
                                      stderr=subprocess.PIPE)
    stdout_wind, stderr_wind = contrarun_wind.communicate()

    names = [os.path.basename(i.split()[4]).split(".bpseq")[0] for i in stdout_unwind.splitlines()]
    winded = [float(i.split()[-1]) for i in stdout_wind.splitlines()]
    unwinded = [float(i.split()[-1]) for i in stdout_unwind.splitlines()]

    for mrid, unw, w in izip(names, unwinded, winded):
        for mirna in mirnas_dict[mrid]:
            ids = mrid.split(",")
            outfile.write("%s,%s,%s,%s\t%f\n" % (ids[0], mirna, ids[1], ids[2], unw - w))
    for bc in bad_coords:
        for mirna in mirnas_dict[bc]:
            ids = bc.split(",")
            outfile.write("%s,%s,%s,%s\t%s\n" % (ids[0], mirna, ids[1], ids[2], "NA"))

    outfile.close()
    shutil.rmtree(dirpath, ignore_errors=True)

if __name__ == '__main__':
    arguments, args = parser.parse_args()
    verbose = arguments.verbose
    main(arguments, verbose)
