#!/usr/bin/env python
import optparse
import sys
import os
import dendropy
import subprocess
import re
import itertools
import traceback
import cPickle as cpickle
import pandas as pd
import itertools as it
from numpy import mean, median, abs
from os.path import abspath
from Bio import SeqIO
from sys import exit

# parse arguments and store them in the variables
parser = optparse.OptionParser()
parser.add_option('-v', '--verbose', dest='verbose', action="store_true",
                     default=False)
parser.add_option('--onlyMIRZA', dest='onlymirza', action="store",
                     default='no', help='calculate only MIRZA score for given coordinates')
parser.add_option('--seq', default='seqs.fa', dest='seq',
                     action='store', help='fasta with mRNA sequences')
parser.add_option('--out', '-o', default='output.tab', dest='out',
                     action='store', help='output table')
parser.add_option('--coords', default='coords.tab', dest='coords',
                     action='store', help='file with target sites positions, miRNA and target gene ID')
parser.add_option('--motifs', default='', dest='motifs',
                     action='store', help='fasta file with miRNA sequences')
parser.add_option('--tree', default='', dest='tree',
                     action='store', help='phylogenetic tree of the species used in alignment file')
parser.add_option('--mln-dir', default='', dest='mln_dir',
                     action='store', help='Directory with multiple alignment files')
parser.add_option('--threshold', type=float, dest='thr', default=20.0,
                     action='store', help='Threshold for MIRZA score')
parser.add_option('--contextLen', type=str, dest='contextLen', default="50",
                     action='store', help='length of the context sequence serounding binding site')
parser.add_option('--mirzabin', dest='mirzabin', action="store",
                     default="", help="Path to the MIRZA binary")
parser.add_option('--reforg', dest='reforg', action="store",
                     default="hg19", help="Reference organism to which alignments are performed: default: hg19")



def main(arguments, verbose):
    #
    # Read coords file into the table: [geneID, miR name, begining position,
    # end position]
    #
    coords = []
    if verbose:
        print "Reading sequences from coordinate file %s" % (arguments.coords)
    try:
        corfile = open(arguments.coords, 'r')
    except IOError:
        raise IOError('Cannot read from coordinate file %s' % (arguments.coords))
    for coord in corfile:
        coord = coord.rstrip().split()
        try:
            coords.append([coord[0], coord[1], int(coord[2]), int(coord[3]), coord[4]])
        except IndexError:
            sys.stderr.write("IndexError in coordinate file, line: %s" % (' '.join(coord)))
        except ValueError, e:
            print str(e)

    #
    # Read mRNA sequences into hash table: {id:sequence}
    #
    if verbose:
        print "Reading sequences from mRNA file %s" % (arguments.seq)
    try:
        seq_obj = open(arguments.seq, 'Ur')
        mRNAseqs = {}
        for seq in SeqIO.parse(seq_obj, 'fasta'):
            mRNAseqs[str(seq.id)] = str(seq.seq)
    except IOError:
        raise IOError('Cannot read from mRNA file %s' % (arguments.seq))

    #
    # In the same way read the miRNA file
    #
    if verbose:
        print "Reading sequences from miRNA file %s" % (arguments.motifs)
    try:
        miRNAseq_obj = open(arguments.motifs, 'Ur')
        miRNAseqs = {}
        for seq in SeqIO.parse(miRNAseq_obj, 'fasta'):
            miRNAseqs[str(seq.id)] = str(seq.seq)
    except IOError:
        raise IOError('Cannot read from miRNA file %s' % (arguments.motifs))


    if arguments.onlymirza != 'yes':
        try:
            phylo_tree = dendropy.Tree()
            phylo_tree.read_from_path(arguments.tree, "newick")
        except IOError:
            raise IOError("Cannot open phylogenetic tree %s" % arguments.tree)

        # make miRNAs for all species
        # get the names of the species in the tree
        species = [leaf.taxon.label for leaf in phylo_tree.leaf_iter()]
        mirhomologues = pd.DataFrame({sp: {mirid: miRNAseqs[mirid][:21]
                                           for mirid in miRNAseqs.keys()}
                                      for sp in species}).transpose()

        # Read multiple alignment file (it is special format - see the exemplary file)
        multiple_alignment_dict = {}
        for coord in coords:
            try:
                handle = open(os.path.join(arguments.mln_dir, coord[0]), 'rb')
                multiple_alignment_dict[coord[0]] = cpickle.load(handle)
            except Exception, e:
                sys.stderr.write("No alignment for %s, going on without it.\n" % coord[0])
                sys.stderr.write(str(e) + "\n")
    #
    # Open output file and write first lines
    #
    try:
        outfile = open(arguments.out, 'w')
    except IOError:
        print "Connot open output file %s" % (arguments.out)

    #
    # write the header to the file
    #
    outfile.write('#siteID\tMIRZAscore\thybrid\tMIRZABranchLengthScoreFill\ttype\tprecise_type\n')

    if verbose is True:
        print "Calculating MIRZA target quality... "


    sys.stderr.write("Collecting sequences\n")
    mRNA_sequences = [cor[-1] for cor in coords]
    mRNA_ids = ["%s,%s,%s" % (cor[0], cor[2], cor[3]) for cor in coords]
    assert len(set([cor[1] for cor in coords])) == 1
    miRNAseq = miRNAseqs[list(set([cor[1] for cor in coords]))[0]][:21]
    miRNAid = list(set([cor[1] for cor in coords]))[0]

    sys.stderr.write("Running MIRZA\n")
    results = calculate_mirza(mRNA_sequences, mRNA_ids, miRNAseq, miRNAid)
    sys.stderr.write("Collecting results\n")
    for key, group in it.groupby(results.splitlines(), lambda x: x == ""):
        if not key:
            proper_group = False
            for line in group:
                if line.startswith(">"):
                    mRNAid = line.split()[0][1:].split(",")[0]
                    beg = line.split()[0][1:].split(",")[1]
                    end = line.split()[0][1:].split(",")[2]
                    score = float(line.split()[-1])
                    proper_group = True
                elif line.startswith("miRNA"):
                    mirhyb = line.split("\t")[1].split(" ")[0]
                elif line.startswith("A L"):
                    hyb = line.split("\t")[1].rstrip()
                elif line.startswith("mRNA"):
                    mrhyb = line.split("\t")[1].split(" ")[0]
            if proper_group:
                hybrids = [mirhyb, hyb, mrhyb]
                mirseq, hybseq, mrhybseq, mrpos = getHybridVector(hybrids)
                canonical, type_of_site = is_canonical([mirseq, hybseq, mrhybseq])
                if arguments.onlymirza != 'yes':
                    try:
                        mln_frag = multiple_alignment_dict[mRNAid]
                        qd = getConservacy(phylotree=phylo_tree,
                                           mrna_frag=mrhyb.replace("-", "")[::-1],
                                           mrnaid=mRNAid,
                                           mirna=mirhomologues,
                                           mirname=miRNAid,
                                           mln_dict=mln_frag,
                                           ref_org=arguments.reforg,
                                           threshold=arguments.thr,
                                           seqlen=int(arguments.contextLen))
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
                outtext = '%s,%s,%s,%s\t%f\t%s\t%s\t%s\t%s\n' % (mRNAid,
                                                             miRNAid,
                                                             beg,
                                                             end,
                                                             score,
                                                             ":".join(hybrids),
                                                             qd,
                                                             "canonical" if canonical else "non-canonical",
                                                             type_of_site)
                outfile.write(outtext)

    outfile.close()
    clean()


def calculate_mirza(mRNAseqs, mRNAids, miRNAseq, miRNAid, update='noupdate'):
    """Get the mirza results

    Args:
        mRNAseqs (@todo): @todo
        mRNAids (@todo): @todo
        miRNAseq (@todo): @todo
        miRNAid (@todo): @todo

    Returns: @todo

    """
    MIRZAbin = arguments.mirzabin
    assert len(mRNAseqs) == len(mRNAids)

    try:
        mrnainput = open(arguments.coords + 'mirza_mrna_input' + '.fa', 'w')
        mirnainput = open(arguments.coords + 'mirza_mirna_input' + '.fa', 'w')
        expressions = open(arguments.coords + 'mirza_mirna_expressions' + '.fa', 'w')
    except IOError:
        raise('Cannot write files for MIRZA input')

    # get the paths for the files
    mrna_path = abspath(mrnainput.name)
    mirna_path = abspath(mirnainput.name)
    expr_path = abspath(expressions.name)

    #
    # write sequences to files
    #
    try:
        for key, value in zip(mRNAids, mRNAseqs):
            if len(value) != 50:
                sys.stderr.write("Skipped %s because wrong length\n" % key)
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



def run_mirza_with_multiple_mrnas(mRNAseq, mRNAID, miRNAseq, miRNAname, conlen, update='noupdate'):
    MIRZAbin = arguments.mirzabin
    #
    # make sure that mirna is 21 in length
    # open all files and write in fasta format
    try:
        mrnainput = open(arguments.coords + 'mirza_mrna_input' + '.fa', 'w')
        mirnainput = open(arguments.coords + 'mirza_mirna_input' + '.fa', 'w')
        expressions = open(arguments.coords + 'mirza_mirna_expressions' + '.fa', 'w')
    except IOError:
        print 'Cannot write files for MIRZA input'
        exit()

    # get the paths for the files
    mrna_path = abspath(mrnainput.name)
    mirna_path = abspath(mirnainput.name)
    expr_path = abspath(expressions.name)

    # write sequences to these files
    try:
        for key, value in mRNAseq.iteritems():
            mrnainput.write('>%s\n%s\n' % (key, value))
        mirnainput.write('>%s\n%s\n' % (miRNAname.replace('-', ''), miRNAseq))
        expressions.write('%s\t1' % (miRNAname))
        mrnainput.close()
        mirnainput.close()
        expressions.close()
    except IOError:
        raise IOError('Cannot write into the MIRZA input files')

    # Run MIRZA with these files
    mirza_command = MIRZAbin + ' ' + expr_path + ' ' + mrna_path + ' ' + mirna_path + ' '\
        + str(len(mRNAseq.values()[0])) + ' ' + update
    mirzarun = subprocess.Popen(mirza_command, shell=True,
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = mirzarun.communicate()
    out_dict = {}
    for line in stdout.split("\n"):
        if line.startswith(">"):
            l = line.split()
            try:
                out_dict[l[0][1:]] = float(l[-1])
            except TypeError, e:
                raise e
    return out_dict


def slideWindows(seq, windowsize, slidesize):
    for i in range(0, len(seq), slidesize):
        # print i,i+windowsize
        if len(seq[i:i + windowsize]) == windowsize:
            yield seq[i:i + windowsize], "-(" + str(i + 1) + "," + str(i + windowsize) + ")"
        else:
            yield seq[-windowsize:], "(" + str(len(seq) - windowsize) + "," + str(len(seq)) + ")"
            break


def clean():
    try:
        os.unlink(arguments.coords + 'mirza_mrna_input' + '.fa')
        os.unlink(arguments.coords + 'mirza_mirna_input' + '.fa')
        os.unlink(arguments.coords + 'mirza_mirna_expressions' + '.fa')
    except:
        pass


def getHybridVector(hyb):
    mirhyb = hyb[0]
    hybrid = hyb[1]
    bpos = re.search("-*[A,C,T,G,U]", hyb[0].upper()).end() - 1
    epos = re.search("-*[A,C,T,G,U]", hyb[0][::-1].upper()).end() - 1
    mrseq = hyb[2][bpos: len(hyb[0]) - epos]
    mirseq = hyb[0][bpos: len(hyb[0]) - epos]
    hybseq = hyb[1][bpos: len(hyb[0]) - epos]
    hybpos = [hybrid.find("|"), hybrid.rfind("|")]
    try:
        assert len(mirhyb) == len(hybrid)
        assert len(mirhyb.replace("-", "")) == 21
    except AssertionError, e:
        print mirseq
        print hybseq
        print mrseq
        raise e
    return mirseq, hybseq, mrseq, hybpos


def extend_both_ends(mrnaseq, span, total_length=50):
    li = span[0]
    ui = span[1]
    length = ui - li
    if length > total_length:
        return -1
    elif length == total_length:
        return mrnaseq[li:ui]
    else:
        dif = total_length - length
        quot = dif // 2  # this is explicit integer division
        l_ext = li - quot
        u_ext = ui + (dif - quot)
        if (l_ext < 0) or (u_ext > len(mrnaseq) - 1):
            return "NA"
        else:
            return mrnaseq[l_ext:u_ext]


def isa_group_separator(line):
    return line == '\n'


def read_mln_to_dict(file, geneslist):
    data = {}
    with open(file) as f:
        for key, block in itertools.groupby(f, isa_group_separator):
            # print key, list(group)
            group = list(block)
            if not key:
                if group[0].rstrip().split(" ")[0][1:].split("|")[0] in geneslist:
                    print group[0]
                    tmp = {}
                    for i in range(0, len(group) - 2, 3):
                        tmp[group[i].rstrip().split(" ")[-1]
                            ] = [group[i + 1].rstrip(), group[i + 2].rstrip()]
                    data[group[i].rstrip().split(" ")[0]
                         [1:].split("|")[0]] = tmp
    return data


def extend_seq(mrnaseq, mrna_frag, total_length=50):
    # Prepare sequences with no gaps
    mrnaseq_nogap = mrnaseq.replace("-", "")
    mrna_frag_nogap = mrna_frag.replace("-", "")
    if len(mrna_frag_nogap) > total_length:
        print "mrnaseq_nogap: ", mrnaseq_nogap
        print "mrna_frag_nogap: ", mrna_frag_nogap
        print "mrnaseq: ", mrnaseq
        print "mrna_frag: ", mrna_frag
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
    # Prepare sequences with no gaps
    mrnaseq_nogap = mrnaseq.replace("-", "")
    mrna_frag_nogap = mrna_frag.replace("-", "")
    if len(mrna_frag_nogap) < total_length:
        print mrna_frag_nogap
        print mrnaseq
        print mrna_frag
        print mrna_frag_target
        raise Exception(
            "Check your sequences maybe you should extend, not shrink them")
        sys.exit()
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


def scale(x, Min, Max):
    return (x - Min) / (Max - Min)


def getConservacy(phylotree, mrna_frag, mrnaid, mirna, mirname, mln_dict,
        seqlen, ref_org="hg19", threshold=50.0):
    # calculate total branch lengths of the tree
    tree = dendropy.Tree(phylotree)
    min_score = dendropy.Tree(phylotree)
    min_score.retain_taxa_with_labels([ref_org])
    total_branch_length = tree.length()
    min_s = min_score.length() / total_branch_length
    reg_pattern = re.compile("-*".join(list(mrna_frag)))
    # get mirna fragments for all organisms
    mrna_len = seqlen
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
    # for each found sequence calculate MIRZAscore
    results_dict = {}
    mrnas_dict = {}
    for horg in queries_dict.keys():
        if queries_dict[horg] != "NA" and mirna[mirname][horg] != "-" * 21:
            # here I should perform checking if miRNA exists in given species
            # score, ht, stdout = runMIRZA(
                # queries_dict[horg], mrnaid, mirna[mirname][horg], mirname, conlen=seqlen)
            mrnas_dict[horg] = queries_dict[horg]
        else:
            results_dict[horg] = "NA"

    if len(mrnas_dict.keys()) > 0:
        results_to_join_dict = run_mirza_with_multiple_mrnas(mrnas_dict,
                                                             mrnaid,
                                                             mirna[mirname][horg],
                                                             mirname,
                                                             conlen=seqlen)
    else:
        results_to_join_dict = {}
    results_dict = dict(results_to_join_dict.items() + results_dict.items())
    # get the list of names where the mirza score is above threshold
    names = [i for i in results_dict.keys() if (
        (results_dict[i] >= threshold) and results_dict[i] != "NA")]
    tree.retain_taxa_with_labels(names + [ref_org])
    # print tree.as_ascii_plot()
    branch_length = tree.length()
    return scale(branch_length / float(total_branch_length), min_s, 1)


def get_indices(desired_length, beg, end):
    """Calculate upper and lower index based on what size we want and what are the coordinates"""
    lower_index = beg - ((desired_length - (end - beg))//2 + (desired_length - (end - beg))%2)
    upper_index = end + ((desired_length - (end - beg))//2)
    return lower_index, upper_index


def is_canonical(hybrids):
    """
    is_cannonical checks if seed is in the mRNA fragment provided and if this
    is the case it marks it as cannonical

    Args:
        hybrids:
            hybrid sequences

    Returns:
        True if the target site is canonical

    Raises:
        pass
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



if __name__ == '__main__':
    arguments, args = parser.parse_args()
    verbose = arguments.verbose
    if arguments.onlymirza not in ["yes", "no"]:
        print "onlyMIRZA option has to be 'yes' or 'no'"
        print parser.print_help()
        sys.exit()
    main(arguments, verbose)
