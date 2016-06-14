#!/usr/bin/env python
"""
Take mirza output and arrange it in proper way
"""

__date_ = "2014-09-16"
__author__ = "Rafal Gumienny"
__email__ = "r.gumienny@unibas.ch"
__license__ = "GPL"

# imports
import re
import sys
import time
import errno
import itertools
from Bio import SeqIO
from argparse import ArgumentParser, RawTextHelpFormatter

parser = ArgumentParser(description=__doc__, formatter_class=RawTextHelpFormatter)
parser.add_argument("-v",
                    "--verbose",
                    dest="verbose",
                    action="store_true",
                    default=False,
                    help="Be loud!")
parser.add_argument("--output",
                    dest="output",
                    help="Output file in Tab format.")
parser.add_argument("--seqs",
                    dest="seqs",
                    required=True,
                    help="UTR sequences in fasta format")
parser.add_argument("--threshold",
                    dest="threshold",
                    type=float,
                    required=True,
                    help="Threshold for the score")
parser.add_argument("--context",
                    dest="context",
                    type=int,
                    default=50,
                    help="Context for sequence to print, defaults to 50")


# redefine a functions for writing to stdout and stderr to save some writting
syserr = sys.stderr.write
sysout = sys.stdout.write


def main(options):
    """Main logic of the script"""
    with open(options.seqs) as f:
        seqs = {}
        for rec in SeqIO.parse(f, 'fasta'):
            seqs[str(rec.id)] = str(rec.seq)

    results = {}
    try:
        for key, group in itertools.groupby(sys.stdin, lambda x: x == "\n"):
            if not key:
                proper_group = False
                for line in group:
                    if line.startswith(">"):
                        mRNAid, begend = line.split()[0][1:].split(":")
                        miRNAid = line.split()[4][1:]
                        beg = int(begend.split(",")[0][1:])
                        end = int(begend.split(",")[1][:-1])
                        score = float(line.split()[-1])
                        if score < options.threshold:
                            break
                        proper_group = True
                    elif line.startswith("miRNA"):
                        mirhyb = line.split("\t")[1].split(" ")[0]
                    elif line.startswith("A L"):
                        hyb = line.split("\t")[1].rstrip()
                    elif line.startswith("mRNA"):
                        mrhyb = line.split("\t")[1].split(" ")[0]
                if proper_group:
                    if score < options.threshold:
                        continue
                    hybrids = [mirhyb, hyb, mrhyb]
                    mirseq, hybseq, mrhybseq, mrpos = get_hybrid_vector(hybrids)
                    canonical, type_of_site = is_canonical([mirseq, hybseq, mrhybseq])
                    # print mRNAid, beg, end, score
                    # print mirseq, hybseq, mrhybseq, mrpos
                    # for h in hybrids:
                    #     print h
                    # print
                    # for h in hybrids:
                    #     print h[::-1]
                    if hybseq[0] == "|":
                        truebeg = end - mrpos[0] - 7
                        trueend = end - mrpos[0]
                    else:
                        truebeg = end - mrpos[0] - 7
                        trueend = end - mrpos[0]
                    # print "Actual position: %i - %i" % (truebeg, trueend)
                    # print seqs[mRNAid][beg - 1: end], len(seqs[mRNAid][beg - 1: end]), seqs[mRNAid][truebeg: trueend]
                    context_beg, context_end = get_indices(options.context, truebeg, trueend)
                    sequence = seqs[mRNAid][context_beg: context_end]
                    if len(sequence) < options.context:
                        continue
                    dkey = "%s\t%s\t%i\t%i" % (mRNAid, miRNAid, truebeg, trueend)
                    dval = (score, "canonical" if canonical else "noncanonical", type_of_site, sequence)
                    if dkey not in results:
                        results[dkey] = dval
                    else:
                        if results[dkey][0] < dval[0]:
                            results[dkey] = dval
    except IOError as e:
        if e.errno == errno.EPIPE:
            pass

    with open(options.output, 'w') as outfile:
        for key, val in results.iteritems():
            outtext = "%s\t%s\n" % (key, "\t".join(str(i) for i in val))
            outfile.write(outtext)

def get_hybrid_vector(hyb):
    """Get the hybrid as defined by miRNA
    binding

    Args:
        hyb (list): hybrids from MIRZA

    Returns: tuple
             ( miRNA sequence (5'-> 3'),
               hybrid shape,
               mRNA sequence (3'-> 5'),
               hybrid position )

    """
    bpos = re.search("-*[A,C,T,G,U]", hyb[0].upper()).end() - 1
    epos = re.search("-*[A,C,T,G,U]", hyb[0][::-1].upper()).end() - 1
    mrseq = hyb[2][bpos: len(hyb[0]) - epos]
    mirseq = hyb[0][bpos: len(hyb[0]) - epos]
    hybseq = hyb[1][bpos: len(hyb[0]) - epos]
    hybpos = [hyb[1].find("|"), hyb[1].rfind("|")]
    try:
        assert len(hyb[0]) == len(hyb[1])
        assert len(hyb[0].replace("-", "")) == 21
    except AssertionError, e:
        syserr(mirseq + "\n")
        syserr(hybseq + "\n")
        syserr(mrseq + "\n")
        # import pdb; pdb.set_trace()
    return mirseq, hybseq, mrseq, hybpos


def is_canonical(hybrids):
    """
    is_cannonical checks if seed is in the mRNA fragment provided and if this
    is the case it marks it as cannonical

    Args:
        hybrids (list): hybrid sequences

    Returns:
        True if the target site is canonical
    """
    mrhyb = hybrids[2].upper().replace("U", "T")
    mirhyb = hybrids[0].upper().replace("U", "T")
    hybrid = hybrids[1]
    """
    2-8
    """
    if hybrid[1:8] == "|||||||":
        guwoble = False
        for mirnuc, mrnuc in zip(mirhyb[1:8], mrhyb[1:8]):
            if (mirnuc == 'G' and mrnuc == 'T') or (mirnuc == 'T' and mrnuc == 'G'):
                guwoble = True
        if guwoble:
            return False, "2-8-Gwoble"
        else:
            return True, "2-8"
    elif (hybrid[1:7] == "||||||" and mrhyb[0] == 'A'):
        guwoble = False
        for mirnuc, mrnuc in zip(mirhyb[1:7], mrhyb[1:7]):
            if (mirnuc == 'G' and mrnuc == 'T') or (mirnuc == 'T' and mrnuc == 'G'):
                guwoble = True
        if guwoble:
            return False, "2-7-A-Gwoble"
        else:
            return True,  "2-7-A"
    else:
        if hybrid[0:7] == "|||||||":
            return False, "1-7-ElMMo"
        elif hybrid[1:7] == "||||||":
            return False, "6-mer"
        if "v" in hybrid[0:8]:
            return False, "mRNAbulge"
        elif "^" in hybrid[0:8]:
            return False, "miRNAbulge"
        elif "O" in hybrid[0:8]:
            return False, "symmetric_loop"
        else:
            return False, "unknown"

def get_indices(desired_length, beg, end):
    """Calculate upper and lower index based on what size we want and what are the coordinates"""
    lower_index = beg - ((desired_length - (end - beg))//2 + (desired_length - (end - beg))%2)
    upper_index = end + ((desired_length - (end - beg))//2)
    return lower_index, upper_index

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
            syserr("############## Started script on %s ##############\n" % start_date)
        main(options)
        if options.verbose:
            syserr("### Successfully finished in %i seconds, on %s ###\n" % (time.time() - start_time, time.strftime("%d-%m-%Y at %H:%M:%S")))
    except KeyboardInterrupt:
        syserr("Interrupted by user after %i seconds!\n" % (time.time() - start_time))
        sys.exit(-1)
