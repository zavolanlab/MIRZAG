#!/usr/bin/env python
"""
Calculate MIRZA and MIRZA BLS for the coordinates provided
"""

__date_ = "2014-05-29"
__author__ = "Rafal Gumienny"
__email__ = "r.gumienny@unibas.ch"
__license__ = "GPL"

# imports
import os
import re
import sys
import time
import gzip
import dendropy
import itertools
import traceback
import subprocess
import cPickle as cpickle
import pandas as pd
from numpy import mean, median, abs
from Bio import SeqIO
from argparse import ArgumentParser

parser = ArgumentParser(description=__doc__)
parser.add_argument("-v",
                    "--verbose",
                    dest="verbose",
                    action="store_true",
                    default=False,
                    help="Be loud!")
parser.add_argument('--onlyMIRZA',
                    dest='onlymirza',
                    choices=('yes', 'no'),
                    default='no',
                    help='Calculate only MIRZA score for given coordinates')
parser.add_argument('--seq',
                    default='seqs.fa',
                    dest='seq',
                    action='store',
                    help='Fasta with mRNA sequences')
parser.add_argument('--out',
                    '-o',
                    default='output.tab',
                    dest='out',
                    action='store',
                    help='Output table')
parser.add_argument('--coords',
                    default='coords.tab',
                    dest='coords',
                    action='store',
                    help='File with target sites positions, miRNA and target gene ID')
parser.add_argument('--motifs',
                    default='',
                    dest='motifs',
                    action='store',
                    help='Fasta file with miRNA sequences')
parser.add_argument('--tree',
                    default='',
                    dest='tree',
                    action='store',
                    help='Phylogenetic tree of the species used in alignment file')
parser.add_argument('--mln-dir',
                    default='',
                    dest='mln_dir',
                    action='store',
                    help='Directory with multiple alignment files')
parser.add_argument('--threshold',
                    type=float,
                    dest='thr',
                    default=20.0,
                    action='store',
                    help='Threshold for MIRZA score')
parser.add_argument('--contextLen',
                    type=int,
                    dest='contextLen',
                    default=50,
                    action='store',
                    help='Length of the context sequence serounding binding site')
parser.add_argument('--mirzabin',
                    dest='mirzabin',
                    action="store",
                    default="",
                    help="Path to the MIRZA binary")
parser.add_argument('--reforg',
                    dest='reforg',
                    action="store",
                    default="hg19",
                    help="Reference organism to which alignments are performed: default: hg19")

try:
    options = parser.parse_args()
except Exception, e:
    parser.print_help()

# redefine a functions for writing to stdout and stderr to save some writting
syserr = sys.stderr.write
sysout = sys.stdout.write


def main():
    """Main logic of the script"""
    if options.verbose:
        syserr("Reading coordinate file\n")
    coords = read_coordinates(options.coords, True)

    if options.verbose:
        syserr("Reading mRNA sequences\n")
    mRNAseqs = read_fasta_to_dict(options.seq)

    if options.verbose:
        syserr("Reading miRNA sequences\n")
    miRNAseqs = read_fasta_to_dict(options.motifs)

    if options.onlymirza != 'yes':
        if options.verbose:
            syserr("Preparing alignments and phylogenetic tree\n")

        phylo_tree = read_phylogenetic_tree(options.tree)
        multiple_alignment_dict = read_multiple_alignments(phylo_tree,
                                                           options.mln_dir,
                                                           coords)
        mirhomologues = make_homologues_mirnas(phylo_tree, miRNAseqs)

    with gzip.open(options.out, 'wb') as outfile:
        if options.verbose:
            syserr("Collecting sequences\n")
        mRNA_sequences = [cor[-1] for cor in coords]
        mRNA_ids = ["%s,%s,%s" % (cor[0], cor[2], cor[3]) for cor in coords]
        number_of_coords = len(set([cor[1] for cor in coords])) == 1
        if number_of_coords > 1:
            raise Exception("More than mirna in coordinate file")
        if number_of_coords == 0:
            syserr("There is no coordinates. Exit.")
            sys.exit()

        miRNAseq = miRNAseqs[list(set([cor[1] for cor in coords]))[0]][:21]
        miRNAid = list(set([cor[1] for cor in coords]))[0]

        if options.verbose:
            syserr("Running MIRZA\n")
        results = calculate_mirza(mRNA_sequences, mRNA_ids, miRNAseq, miRNAid)

        if options.verbose:
            syserr("Collecting results\n")
        for key, group in itertools.groupby(results.splitlines(), lambda x: x == ""):
            if not key:
                proper_group = False
                for line in group:
                    if line.startswith(">"):
                        mRNAid = line.split()[0][1:].split(",")[0]
                        beg = line.split()[0][1:].split(",")[1]
                        end = line.split()[0][1:].split(",")[2]
                        score = float(line.split()[-1])
                        proper_group = True
                    # elif line.startswith("miRNA"):
                    #     mirhyb = line.split("\t")[1].split(" ")[0]
                    # elif line.startswith("A L"):
                    #     hyb = line.split("\t")[1].rstrip()
                    elif line.startswith("mRNA"):
                        mrhyb = line.split("\t")[1].split(" ")[0]
                if proper_group:
                    if len(miRNAseq) < 21:
                        outtext = '%s,%s,%s,%s\t%s\t%s\t%s\t%s\t%s\n' % (mRNAid,
                                                                         miRNAid,
                                                                         beg,
                                                                         end,
                                                                         "NA",
                                                                         "NA",
                                                                         "NA",
                                                                         "NA",
                                                                         "NA")
                        outfile.write(outtext)
                        continue

                    # hybrids = [mirhyb, hyb, mrhyb]
                    # mirseq, hybseq, mrhybseq, mrpos = get_hybrid_vector(hybrids)
                    # canonical, type_of_site = is_canonical([mirseq, hybseq, mrhybseq])
                    if options.onlymirza != 'yes':
                        try:
                            mln_frag = multiple_alignment_dict[mRNAid]
                            qd = calculate_conservation(phylotree=phylo_tree,
                                                        mrna_frag=mrhyb.replace("-", "")[::-1],
                                                        mrnaid=mRNAid,
                                                        mirna=mirhomologues,
                                                        mirname=miRNAid,
                                                        mln_dict=mln_frag,
                                                        ref_org=options.reforg,
                                                        threshold=options.thr,
                                                        mrna_len=options.contextLen)
                            qd = str(qd)
                        except KeyError, e:
                            qd = "NA"
                            sys.stderr.write("KeyError:  " + str(e) + "\n")
                            sys.stderr.write("Trace:  "
                                             + traceback.format_exc()
                                             + "\n")
                            # raise KeyError
                    else:
                        qd = "NA"
                    outtext = '%s,%s,%s,%s\t%f\t%s\n' % (mRNAid,
                                                                     miRNAid,
                                                                     beg,
                                                                     end,
                                                                     score,
                                                                     # ":".join(hybrids),
                                                                     qd)
                                                                     # "canonical" if canonical else "non-canonical",
                                                                     # type_of_site)
                    outfile.write(outtext)
    clean()



def calculate_mirza(mRNAseqs, mRNAids, miRNAseq, miRNAid, update='noupdate',
        context_len=50):
    """Get the mirza results

    Args:
        mRNAseqs (@todo): @todo
        mRNAids (@todo): @todo
        miRNAseq (@todo): @todo
        miRNAid (@todo): @todo

    Kwargs:
        update (str): update prior probability? defaults to noupdate

    Returns: @todo

    """
    MIRZAbin = options.mirzabin
    assert len(mRNAseqs) == len(mRNAids)

    try:
        mrnainput = open(options.coords + 'mirza_mrna_input' + '.fa', 'w')
        mirnainput = open(options.coords + 'mirza_mirna_input' + '.fa', 'w')
        expressions = open(options.coords + 'mirza_mirna_expressions' + '.fa', 'w')
    except IOError:
        raise('Cannot write files for MIRZA input')

    # get the paths for the files
    mrna_path = os.path.abspath(mrnainput.name)
    mirna_path = os.path.abspath(mirnainput.name)
    expr_path = os.path.abspath(expressions.name)

    #
    # write sequences to files
    #
    try:
        for key, value in zip(mRNAids, mRNAseqs):
            if len(value) != context_len:
                sys.stderr.write("Skipped %s because of wrong length\n" % key)
                continue
            mrnainput.write('>%s\n%s\n' % (key, value))
        mirnainput.write('>%s\n%s\n' % (miRNAid.replace('-', ''), miRNAseq))
        expressions.write('%s\t1' % (miRNAid))
        mrnainput.close()
        mirnainput.close()
        expressions.close()
    except IOError:
        raise IOError('Cannot write into the MIRZA input files')
    #
    # Run MIRZA with these files
    #
    mirza_command = MIRZAbin + ' ' + expr_path + ' ' + mrna_path + ' ' + mirna_path + ' '\
        + str(len(mRNAseqs[0])) + ' ' + update
    mirzarun = subprocess.Popen(mirza_command, shell=True,
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = mirzarun.communicate()
    return stdout


def calculate_mirza_for_mln(mRNAseqs, mRNAids, miRNAseq, miRNAid,
        update='noupdate', context_len=50):
    """Get the mirza results

    Args:
        mRNAseqs (@todo): @todo
        mRNAids (@todo): @todo
        miRNAseq (@todo): @todo
        miRNAid (@todo): @todo

    Kwargs:
        update (str): update prior probability? defaults to noupdate

    Returns: @todo

    """
    #
    # Get MIRZA results
    #
    stdout = calculate_mirza(mRNAseqs,
                             mRNAids,
                             miRNAseq,
                             miRNAid,
                             update,
                             context_len)

    #
    # Parse results into dictionary
    out_dict = {}
    for line in stdout.split("\n"):
        if line.startswith(">"):
            l = line.split()
            try:
                out_dict[l[0][1:]] = float(l[-1])
            except TypeError, e:
                raise e
    return out_dict


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


def extend_seq(mrnaseq, mrna_frag, total_length=50):
    """Extend the sequence from both ends

    Args:
        mrnaseq (str): sequence (with gaps) to extend from
        mrna_frag (str): sequence (with gaps) for extension

    Kwargs:
        total_length (int): length to extend to

    Returns: @todo

    """
    #
    # Prepare sequences with no gaps
    #
    mrnaseq_nogap = mrnaseq.replace("-", "")
    mrna_frag_nogap = mrna_frag.replace("-", "")
    #
    # check if the sequence is shorter
    #
    if len(mrna_frag_nogap) > total_length:
        syserr("mrnaseq_nogap: ", mrnaseq_nogap)
        syserr("mrna_frag_nogap: ", mrna_frag_nogap)
        syserr("mrnaseq: ", mrnaseq)
        syserr("mrna_frag: ", mrna_frag)
        raise Exception(
            "Check your sequences maybe you should shrink, not extend them")
    span = re.search(mrna_frag_nogap, mrnaseq_nogap).span()

    # Decide which type of extension to do
    gap_pos_mean = mean([i for i, x in enumerate(mrna_frag) if x == "-"])
    list_median = median([i for i in range(len(mrna_frag))])

    # this ratio gives us relative position of the gaps
    ratio = gap_pos_mean / list_median

    # Based on the ratio do the extension of the sequence
    if ratio > 0.5 and ratio < 1.5:  # extend both sides
        li = span[0]
        ui = span[1]
        length = ui - li
        if length > total_length:
            return -1
        elif length == total_length:
            return mrnaseq_nogap[li:ui]
        else:
            dif = total_length - length
            quot = dif // 2  # this is explicit integer division
            l_ext = li - quot  # TODO check if they are not lower than 0
            u_ext = ui + (dif - quot)
            if (l_ext < 0) or (u_ext > len(mrnaseq_nogap) - 1):
                return "NA"
            else:
                return mrnaseq_nogap[l_ext:u_ext]
    elif ratio <= 0.5:  # extend left - it means upstream (5'end)
        li = span[0]
        ui = span[1]
        length = ui - li
        dif = total_length - len(mrna_frag_nogap)
        if (li - dif < 0):
            return mrnaseq_nogap[:ui + abs(li - dif)]
        else:
            return mrnaseq_nogap[li - dif:ui]
    elif ratio >= 1.5:  # extend right - it means downstream (3'end)
        li = span[0]
        ui = span[1]
        length = ui - li
        dif = total_length - len(mrna_frag_nogap)
        # if there is noting to extend to the right
        if ui + dif > len(mrnaseq_nogap):
            return mrnaseq_nogap[li - ((ui + dif) - len(mrnaseq_nogap)):]
        else:
            return mrnaseq_nogap[li:ui + dif]


def shrink_seq(mrnaseq, mrna_frag, mrna_frag_target, total_length=50):
    """Extend the sequence from both ends

    Args:
        mrnaseq (str): sequence (with gaps) to extend from
        mrna_frag (str): sequence (with gaps) for extension
        mrna_frag_target (str): sequence of the target mRNA

    Kwargs:
        total_length (int): length to extend to

    Returns: @todo

    """
    # Prepare sequences with no gaps
    mrnaseq_nogap = mrnaseq.replace("-", "")
    mrna_frag_nogap = mrna_frag.replace("-", "")
    if len(mrna_frag_nogap) < total_length:
        syserr(mrna_frag_nogap)
        syserr(mrnaseq)
        syserr(mrna_frag)
        syserr(mrna_frag_target)
        raise Exception(
            "Check your sequences maybe you should extend, not shrink them")
    span = re.search(mrna_frag_nogap, mrnaseq_nogap).span()

    # Decide which type of extension to do
    gap_pos_mean = mean(
        [i for i, x in enumerate(mrna_frag_target) if x == "-"])
    list_median = median([i for i in range(len(mrna_frag_target))])

    # this ratio gives us relative position of the gaps
    ratio = gap_pos_mean / list_median

    # Based on the ratio do the shrinkage of the sequence
    if ratio > 0.5 and ratio < 1.5:  # extend both sides
        li = span[0]
        ui = span[1]
        length = ui - li
        if length < total_length:
            return -1
        elif length == total_length:
            return mrnaseq_nogap[li:ui]
        else:
            dif = abs(total_length - length)
            quot = dif // 2  # this is explicit integer division
            l_ext = li + quot
            u_ext = ui - (dif - quot)
            if (u_ext < 0) or (u_ext > len(mrnaseq_nogap) - 1):
                return "NA"
            else:
                return mrnaseq_nogap[l_ext:u_ext]
    elif ratio <= 0.5:  # trim left - it means upstream (5'end)
        li = span[0]
        ui = span[1]
        length = ui - li
        dif = len(mrna_frag_nogap) - total_length
        return mrnaseq_nogap[li + dif:ui]
    elif ratio >= 1.5:  # extend right - it means downstream (3'end)
        li = span[0]
        ui = span[1]
        length = ui - li
        dif = len(mrna_frag_nogap) - total_length
        return mrnaseq_nogap[li:ui - dif]


def scale(x, minimum, maximum):
    """Scale number from 0 to 1

    Args:
        x (float): number to scale
        minimum (float): minimum for the set
        maximum (float): maximum for the set

    Returns: float

    """
    return (x - minimum) / (maximum - minimum)


def calculate_conservation(phylotree, mrna_frag, mrnaid, mirna, mirname,
                           mln_dict, mrna_len, ref_org='hg19', threshold=50):
    """Calculate conservation

    Args:
        phylotree (@todo): @todo
        mrna_frag (@todo): @todo
        mrnaid (@todo): @todo
        mirna (@todo): @todo
        mirname (@todo): @todo
        mln_dict (@todo): @todo
        mrna_len (@todo): @todo

    Kwargs:
        ref_org (@todo): @todo
        threshold (@todo): @todo

    Returns: float

    """
    #
    # calculate total branch lengths of the tree
    #
    tree = dendropy.Tree(phylotree)
    min_score = dendropy.Tree(phylotree)
    min_score.retain_taxa_with_labels([ref_org])
    total_branch_length = tree.length()
    min_s = min_score.length() / total_branch_length
    #
    # get the dictionary with the homologous fragments
    #
    reg_pattern = re.compile("-*".join(list(mrna_frag)))
    queries_dict = {}
    for org in mln_dict.keys():
        TSeq = mln_dict[org][0]
        QSeq = mln_dict[org][1]
        match = re.search(reg_pattern, TSeq)
        if match != None and len(QSeq.replace("-", "")) >= mrna_len:
            span = match.span()
            tseq = TSeq[span[0]:span[1]]  # those are sequences with gaps
            qseq = QSeq[span[0]:span[1]]
            # Decide if we need to shrink or expand sequences
            if len(qseq.replace("-", "")) < mrna_len:
                myseq = extend_seq(QSeq, qseq, total_length=mrna_len)
                queries_dict[org] = myseq
                try:
                    assert len(myseq) == mrna_len or myseq == "NA"
                except AssertionError:
                    message = "Error in extend_seq"
                    raise Exception(message)
            elif len(qseq.replace("-", "")) > mrna_len:
                myseq = shrink_seq(QSeq, qseq, tseq, total_length=mrna_len)
                queries_dict[org] = myseq
                assert len(myseq) == mrna_len or myseq == "NA"
            else:
                if "-" in qseq:
                    myseq = qseq.replace("-", "")
                    assert len(myseq) == mrna_len or myseq == "NA"
                else:
                    myseq = qseq
                queries_dict[org] = myseq
        else:
            queries_dict[org] = "NA"
            myseq = "NA"
    #
    # for each found sequence calculate MIRZAscore
    #
    results_dict = {}
    mrnas_dict = {}
    for horg in queries_dict.keys():
        if queries_dict[horg] != "NA":  # and mirna[mirname][horg] != "-" * 21:
            mrnas_dict[horg] = queries_dict[horg]
        else:
            results_dict[horg] = "NA"

    if len(mrnas_dict.keys()) > 0:
        results_to_join_dict = calculate_mirza_for_mln(mrnas_dict.values(),
                                                       mrnas_dict.keys(),
                                                       mirna[mirname][horg],
                                                       mirname,
                                                       update='noupdate',
                                                       context_len=mrna_len)
    else:
        results_to_join_dict = {}
    #
    # Join the dictionaries with scores and with NAs
    #
    results_dict = dict(results_to_join_dict.items() + results_dict.items())

    # get the list of names where the mirza score is above threshold
    names = [i for i in results_dict.keys() if
             ((results_dict[i] >= threshold) and results_dict[i] != "NA")]
    tree.retain_taxa_with_labels(names + [ref_org])
    branch_length = tree.length()
    conservation = scale(branch_length / float(total_branch_length), min_s, 1)
    return conservation


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


def read_coordinates(path_to_file, zipped=False):
    """Read coordinate file into list

    Args:
        path_to_file (@todo): @todo

    Returns: list of lists

    """
    coords = []
    if options.verbose:
        syserr("Reading sequences from coordinate file %s\n" % (path_to_file))
    try:
        if zipped:
            corfile = gzip.open(path_to_file, 'rb')
        else:
            corfile = open(path_to_file, 'r')
    except IOError:
        raise IOError('Cannot read from coordinate file %s' % (path_to_file))
    for coord in corfile:
        coord = coord.rstrip().split()
        try:
            coords.append([coord[0], coord[1], int(coord[2]), int(coord[3]), coord[4]])
        except IndexError:
            syserr("IndexError in coordinate file, line: %s\n" % (' '.join(coord)))
        except ValueError, e:
            syserr(str(e) + "\n")

    return coords


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


def read_phylogenetic_tree(path_to_file):
    """Read phylogenetic tree

    Args:
        path_to_file (@todo): @todo

    Returns: @todo

    """
    try:
        phylo_tree = dendropy.Tree()
        phylo_tree.read_from_path(path_to_file, "newick")
    except IOError:
        raise IOError("Cannot open phylogenetic tree %s" % path_to_file)

    return phylo_tree


def read_multiple_alignments(tree, directory_path, coordinates):
    """Read MLN files into dict

    Args:
        tree (@todo): @todo
        directory_path (@todo): @todo
        coordinates (@todo): @todo

    Returns: @todo

    """
    multiple_alignment_dict = {}
    for coord in coordinates:
        try:
            handle = open(os.path.join(directory_path, coord[0]), 'rb')
            multiple_alignment_dict[coord[0]] = cpickle.load(handle)
        except Exception, e:
            syserr("No alignment for %s, going on without it.\n" % coord[0])
            syserr(str(e) + "\n")

    return multiple_alignment_dict


def make_homologues_mirnas(phylogenetic_tree, mirna_seqs):
    """Make DataFrame with miRNAs for different species

    Args:
        phylogenetic_tree (@todo): @todo
        mirna_seqs (@todo): @todo

    Returns: @todo

    """
    species = [leaf.taxon.label for leaf in phylogenetic_tree.leaf_iter()]
    mirhomologues = pd.DataFrame({sp: {mirid: mirna_seqs[mirid][:21]
                                       for mirid in mirna_seqs.keys()}
                                  for sp in species}).transpose()
    return mirhomologues


def clean():
    """Remove files after calculation
    """
    try:
        os.unlink(options.coords + 'mirza_mrna_input' + '.fa')
        os.unlink(options.coords + 'mirza_mirna_input' + '.fa')
        os.unlink(options.coords + 'mirza_mirna_expressions' + '.fa')
    except:
        pass

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
