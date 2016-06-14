#!/usr/bin/env python
"""
Split coordinate file by miRNA
"""

__date_ = "2014-09-15"
__author__ = "Rafal Gumienny"
__email__ = "r.gumienny@unibas.ch"
__license__ = "GPL"

# imports
import sys
import time
from collections import defaultdict
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
                    help="Input file in Tab format.")
parser.add_argument("--prefix",
                    dest="prefix",
                    default="./",
                    help="Prefix to the file, defaults to ./")
parser.add_argument("--suffix",
                    dest="suffix",
                    default=".coords",
                    help="Suffix to the file, defaults to .coords")
parser.add_argument("--column",
                    dest="column",
                    type=int,
                    default=1,
                    help="Keyword in column, defaults to 1")
parser.add_argument("--gzip",
                    dest="gzip",
                    action="store_true",
                    default=False,
                    help="Gzip output file")

try:
    options = parser.parse_args()
except Exception, e:
    parser.print_help()

# redefine a functions for writing to stdout and stderr to save some writting
syserr = sys.stderr.write
sysout = sys.stdout.write

if options.gzip:
    import gzip

# class MyDict(dict):

#     """A dictionary object to store open files"""

#     def __init__(self, factory):
#         dict.__init__(self)
#         self.factory = factory

#     def __missing__(self, key):
#         self[key] = self.factory(key)
#         return self[key]


def main():
    """Main logic of the script"""
    # if options.gzip:
    #     files = MyDict(lambda x: gzip.open(x, "a"))
    # else:
    #     files = MyDict(lambda x: open(x, "a"))

    if options.gzip:
        with open(options.input) as infile:
            for line in infile:
                filename = options.prefix + line.split("\t")[options.column] + options.suffix
                with gzip.open(filename, 'ab') as outfile:
                    outfile.write(line)
    else:
        with open(options.input) as infile:
            for line in infile:
                filename = options.prefix + line.split("\t")[options.column] + options.suffix
                with open(filename, 'a') as outfile:
                    outfile.write(line)

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
