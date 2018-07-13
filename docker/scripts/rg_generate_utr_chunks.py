#!/usr/bin/env python
"""
Take UTRs and generate fragments by sliding window that can be fed into MIRZA
"""

__date__ = "2015-06-15"
__author__ = "Rafal Gumienny"
__email__ = "r.gumienny@unibas.ch"
__license__ = "GPL"

# imports
import os
import sys
import time
from Bio import SeqIO
from contextlib import contextmanager
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
                    default=sys.stdin,
                    help="Input file in fasta format. Defaults to sys.stdin.")
parser.add_argument("--output-dir",
                    dest="output_dir",
                    default="",
                    help="Output directory for split files")
parser.add_argument("--part-size",
                    dest="part_size",
                    type=int,
                    default=40000,
                    help="Number of sequences per part, defaults to 40000")
parser.add_argument("--window-size",
                    dest="window_size",
                    type=int,
                    default=50,
                    help="Length of the window for MIRZA, defaults to 50")
parser.add_argument("--slide-size",
                    dest="slide_size",
                    type=int,
                    default=20,
                    help="Size of the window slide, defaults to 20")


# redefine a functions for writing to stdout and stderr to save some writting
syserr = sys.stderr.write
sysout = sys.stdout.write


def main(options):
    """Main logic of the script"""
    try:
        size_count = 0
        outfile = None
        size_count = 0
        files_count = 0
        with smart_open(options.input) as infile:
            for rec in SeqIO.parse(infile, 'fasta'):
                windows = slide_windows(str(rec.seq),
                                        options.window_size,
                                        options.slide_size)
                for part_mrna, winpos, winpos_tup in windows:
                    if size_count % options.part_size == 0:
                        if outfile is not None:
                            outfile.close()
                        files_count += 1
                        outfile = open(os.path.join(options.output_dir, "part_%i.mrna.fa" % files_count), 'w')
                    if len(part_mrna) == options.window_size:
                        outfile.write(">%s:%s\n%s\n" % (str(rec.id), winpos, part_mrna))
                        size_count += 1
    except Exception, e:
        raise e
    finally:
        if outfile is not None:
            outfile.close()

    if options.verbose:
        syserr("File was split into %i chunks.\n" % files_count)


# this function is also defined in utils but I put it here to avoid
# unnecessary module import that might not be available everywhere as
# it is my own module
@contextmanager
def smart_open(filepath, mode='r'):
    """Open file intelligently depending on the source

    :param filepath: can be both path to file or sys.stdin or sys.stdout
    :param mode: mode can be read "r" or write "w". Defaults to "r"
    :yield: context manager for file handle

    """
    if mode == 'r':
        if filepath is not sys.stdin:
            fh = open(filepath, 'r')
        else:
            fh = filepath
        try:
            yield fh
        except IOError as e:
            if fh is not sys.stdin:
                fh.close()
            elif e.errno == errno.EPIPE:
                pass
        finally:
            if fh is not sys.stdin:
                fh.close()
    elif mode == 'w':
        if filepath is not sys.stdout:
            fh = open(filepath, 'w')
        else:
            fh = filepath
        try:
            yield fh
        finally:
            if fh is not sys.stdout:
                fh.close()
    else:
        raise NoSuchModeException("No mode %s for file" % mode)


def slide_windows(seq, windowsize, slidesize):
    """
    Slide window on string in order to divide it into parts of equal windowsize length
    """
    for i in range(0, len(seq), slidesize):
        if len(seq[i:i+windowsize]) == windowsize:
            # yield a part of the sequence of length windowsize
            yield seq[i:i+windowsize], "("+str(i+1)+","+str(i+windowsize)+")", (i+1, i+windowsize)
        else:
            # if we reach end of the sequence return last fragment of size 50
            yield seq[-windowsize:], "("+str(len(seq)-windowsize)+","+str(len(seq))+")", (len(seq)-windowsize, len(seq))
            break


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
