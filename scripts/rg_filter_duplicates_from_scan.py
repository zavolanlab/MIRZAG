#!/usr/bin/env python
"""
Filter duplicates in coordinates by id, miRNA and sequence
"""

__date_ = "2014-06-15"
__author__ = "Rafal Gumienny"
__email__ = "r.gumienny@unibas.ch"
__license__ = "GPL"

# imports
import sys
import time
import gzip
import pandas as pd
from argparse import ArgumentParser

parser = ArgumentParser(description=__doc__)
parser.add_argument("-v",
                    "--verbose",
                    dest="verbose",
                    action="store_true",
                    default=False,
                    help="Be loud!")
parser.add_argument("--coords",
                    dest="coords",
                    required=True,
                    help="Coordinates")
parser.add_argument("--split-by",
                    dest="split_by",
                    required=True,
                    help="Split id by the string")
parser.add_argument("--index-after-split",
                    dest="index_after_split",
                    type=int,
                    required=True,
                    help="After split take this column as new id, 0 based")
parser.add_argument("--output",
                    dest="output",
                    required=True,
                    help="Output name")


# redefine a functions for writing to stdout and stderr to save some writting
syserr = sys.stderr.write
sysout = sys.stdout.write


def main(options):
    """Main logic of the script"""
    names = ['id', 'mirna', 'beg', 'end','score','type','ptype', 'seq']
    names_to_save = ['id', 'mirna', 'beg', 'end', 'seq']
    df = pd.read_table(options.coords, names=names)
    df['newid'] = [i.split(options.split_by)[options.index_after_split] for i in df.id]
    ndf = df.drop_duplicates(cols=['newid', 'mirna', 'seq'])
    with open(options.output, 'wb') as o:
        ndf[names_to_save].to_csv(o, header=False, index=False, sep='\t')
    if options.verbose:
        syserr("Filtered %s\n" % options.coords)
        syserr(" - number of coordinates before: %i\n" % len(df.index))
        syserr(" - number of coordinates after : %i\n" % len(ndf.index))

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
