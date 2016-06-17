#!/usr/bin/env python
"""
A `CONTRAfold <http://contra.stanford.edu/contrafold/>`_ algorithm is used to calculate accessibility of target site for miRNA.
"""

__date__ = "2014-12-19"
__author__ = "Rafal Gumienny"
__email__ = "r.gumienny@unibas.ch"
__license__ = "GPL"

# imports
import sys
import time
import os
import gzip
import shutil
import optparse
import subprocess
from Bio import SeqIO
from itertools import izip
from argparse import ArgumentParser, RawTextHelpFormatter

parser = ArgumentParser(description=__doc__, formatter_class=RawTextHelpFormatter)
parser.add_argument("-v",
                    "--verbose",
                    dest="verbose",
                    action="store_true",
                    default=False,
                    help="Be loud!")
parser.add_argument("--seq",
                    dest="seq",
                    default="seqs.fa",
                    help="Fasta file with mRNA sequences , defaults to seqs.fa")
parser.add_argument("--out",
                    dest="out",
                    default="output.tab.gz",
                    help="output file, defaults to output.tab.gz")
parser.add_argument('--coords',
                  default='coords.tab',
                  dest='coords',
                  action='store',
                  help='file with target\
                  sites positions, miRNA and target gene ID')
parser.add_argument('--contextLen_L',
                  type=int,
                  dest='contextLen_L',
                  default=0,
                  action='store',
                  help='length of the \
                  context sequence downstream binding site to be unwinded')
parser.add_argument('--contextLen_U',
                  type=int,
                  dest='contextLen_U',
                  default=0,
                  action='store',
                  help='length of the \
                  context sequence upstream binding site to be unwinded')
parser.add_argument('--context',
                  type=int,
                  dest='context',
                  default=50,
                  action='store',
                  help='length of the \
                  context of the seed to be checked')
parser.add_argument('--contrabin',
                     default="contrafold",
                     help="Path to CONTRAfold binary")


# redefine a functions for writing to stdout and stderr to save some writting
syserr = sys.stderr.write
sysout = sys.stdout.write


def main(options):
    """Main logic of the script"""
    # Check if path to CONTRAfold is valid
    if not is_executable(options.contrabin):
        raise Exception("Path to CONTRAfold is invalid (%s)! Please define it with --contrabin option." % options.contrabin)
    # Read coords file into the table: [geneID, miR name, begining position,
    # end position]
    coords = []
    if options.verbose is True:
        print "Reading sequences from coordinate file %s" % (options.coords)
    try:
        corfile = gzip.open(options.coords)
    except IOError:
        raise IOError('Cannot read from coordinate file %s' % (options.coords))

    for c in corfile:
        try:
            c = c.rstrip().split()
            coords.append([c[0], int(c[2]), int(c[3]), c[1]])
        except ValueError:
            raise ValueError("Wrong coordinates: %s" % " ".join(c))


    # Read mRNA sequences into hash table: {id:sequence}
    if options.verbose is True:
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
    # outfile.write('#siteID\tContraScoreTargetSite\n')

    # Iterate through the binding coordinates to calculate their score
    if options.verbose is True:
        print "Calculating Contrafold accessibility..."
    count = 0
    corlen = len(coords)

    # make a directory to store the files
    dirpath = os.path.splitext(os.path.abspath(options.coords))[0]
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
            lowerIndex = options.context - options.contextLen_L
            upperIndex = options.context + options.contextLen_U + coorLen
            mrna_fragment = mrnasequ[lowerix - options.context:upperix + options.context]
            mrnarange = range(1, len(mrna_fragment) + 1)
            if len(mrnasequ) >= upperIndex and  lowerIndex >= 0 \
                    and len(mrna_fragment) == 2 * options.context + coorLen:
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
                            (options.contrabin,
                             os.path.join(dirpath, "*.bpseq"))
    contrarun_unwind = subprocess.Popen(contra_command_unwind,
                                        shell=True,
                                        stdout=subprocess.PIPE,
                                        stderr=subprocess.PIPE)
    stdout_unwind, stderr_unwind = contrarun_unwind.communicate()

    contra_command_wind = "%s predict %s --partition" % \
                          (options.contrabin,
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


def is_executable(program):
    """
    Check if the path/binary provided is valid executable
    """
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return True
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return True

    return False

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
