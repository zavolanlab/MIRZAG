#!/usr/bin/env python
"""
Prepare miRNA fasta file for MIRZA i.e. for each miRNA sequence
check if it is 21 nucleotide long and if not eliminate it. It
also replaces all u or U into T. In the same time it splits miRNAs
into separate files.
"""

__date__ = "2014-12-13"
__author__ = "Rafal Gumienny"
__email__ = "r.gumienny@unibas.ch"
__license__ = "GPL"

# imports
import os
import sys
import time
from Bio import SeqIO
from argparse import ArgumentParser, RawTextHelpFormatter


parser = ArgumentParser(description=__doc__, formatter_class=RawTextHelpFormatter)
parser.add_argument("-v",
                    "--verbose",
                    dest="verbose",
                    action="store_true",
                    default=False,
                    help="Be loud!")
parser.add_argument("--input",
                    dest="input",
                    required=True,
                    help="Input miRNA file in fasta format.")
parser.add_argument("--output-dir",
                    dest="output_dir",
                    default="",
                    help="Directory for split files, defaults to Output")


# redefine a functions for writing to stdout and stderr to save some writting
syserr = sys.stderr.write
sysout = sys.stdout.write


def main(options):
    """Main logic of the script"""
    with open(options.input) as infile:
        for rec in SeqIO.parse(infile, 'fasta'):
            if len(rec.seq) < 21:
                syserr("mi/siRNA %s is shorter than 21 nucleotides. " % (rec.id) + \
                        "It will be removed from the list.\n")
            else:
                with open(os.path.join(options.output_dir, str(rec.id) + ".mirna.fa"), 'w') as outfile:
                    outfile.write(">%s\n%s\n" % (rec.id,
                                                 str(rec.seq)[:21].upper().replace("U", "T")))


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
