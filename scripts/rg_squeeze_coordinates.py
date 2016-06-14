#!/usr/bin/env python
"""

"""

__date_ = "2014-06-19"
__author__ = "Rafal Gumienny"
__email__ = "r.gumienny@unibas.ch"
__license__ = "GPL"

# imports
import sys
import csv
import time
import gzip
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
                    help="Input coordinate file")
parser.add_argument("--output",
                    dest="output",
                    help="output.tab")


# redefine a functions for writing to stdout and stderr to save some writting
syserr = sys.stderr.write
sysout = sys.stdout.write


def main(options):
    """Main logic of the script"""
    data = defaultdict(list)
    with open(options.input) as f:
        for sid, mirna, beg, end, seq in csv.reader(f, delimiter='\t'):
            data[','.join([sid, beg])].append(mirna)

    with gzip.open(options.output, 'wb') as o:
        for key, value in data.iteritems():
            sid, beg = key.split(",")
            beg = int(beg)
            o.write("%s\t%i\t%i\t%s\n" % (sid, beg, beg + 7, ",".join(value)))

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
