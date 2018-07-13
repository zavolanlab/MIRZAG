#!/usr/bin/env python
"""
Collect all sites and calculate probability
"""

__date_ = "2014-07-14"
__author__ = "Rafal Gumienny"
__email__ = "r.gumienny@unibas.ch"
__license__ = "GPL"

# imports
import os
import sys
import csv
import gzip
import time
from collections import defaultdict
from argparse import ArgumentParser

parser = ArgumentParser(description=__doc__)
parser.add_argument("-v",
                    "--verbose",
                    dest="verbose",
                    action="store_true",
                    default=False,
                    help="Be loud!")
parser.add_argument("--input",
                    dest="input",
                    required=True,
                    help="Input file")
parser.add_argument("--output",
                    dest="output",
                    default="output.tab",
                    help="Output name , defaults to output.tab")
parser.add_argument("--threshold",
                    dest="threshold",
                    type=float,
                    default=0.12,
                    help="Threshold for summing, defaults to 0.12")
parser.add_argument("--split-by",
                    dest="split_by",
                    default="|",
                    help="If the header of fasta has multiple annotations eg. transcript_id|entrez_id|wikiname\nsplit it and take only one, defaults to |")
parser.add_argument("--column",
                    dest="column",
                    type=int,
                    default=1,
                    help="0 based column number to take after splitting, defaults to 1")
parser.add_argument("--name",
                    dest="name",
                    default="GeneID",
                    help="Name of the id eg. gene, transcript etc , defaults to GeneID")


# redefine a functions for writing to stdout and stderr to save some writting
syserr = sys.stderr.write
sysout = sys.stdout.write


def main(options):
    """Main logic of the script"""
    data_with_conserved = {}
    data_without_conserved = {}
    if options.verbose:
        syserr("Reading Data:\n")
    with gzip.open(options.input) as infile:
        for row in csv.reader(infile, delimiter='\t'):
            if row[0] != 'ID':
                key = "%s,%s" % (row[0].split(options.split_by)[options.column], row[1])
                if row[10] != 'NaN':
                    if float(row[10]) >= options.threshold:
                        if key not in data_without_conserved:
                            data_without_conserved[key] = 0.0
                        data_without_conserved[key] += float(row[10])
                if row[11] != 'NaN':
                    if float(row[11]) >= options.threshold:
                        if key not in data_with_conserved:
                            data_with_conserved[key] = 0.0
                        data_with_conserved[key] += float(row[11])
    if options.verbose:
        syserr("Writing output\n")
    with gzip.open(options.output, 'wb') as o:
        # o.write("%s\t%s\t%s\t%s\n" % (options.name, 'miRNA', 'Total score without conservation', 'Total score with conservation'))
        for key, value in data_without_conserved.iteritems():
            myid, mirna = key.split(",")
            try:
                with_conservation_value = str(data_with_conserved[key])
            except KeyError:
                with_conservation_value = 'NaN'
            o.write("%s\t%s\t%f\t%s\n" % (myid, mirna, value, with_conservation_value))


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
