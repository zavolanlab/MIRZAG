#!/usr/bin/env python
"""
Merge all results into one features table
"""

__date_ = "2014-06-05"
__author__ = "Rafal Gumienny"
__email__ = "r.gumienny@unibas.ch"
__license__ = "GPL"

# imports
import gc
import sys
import time
import pandas as pd
import numpy as np
import statsmodels.api as sm
from itertools import starmap, izip
from operator import mul
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
                    help="Input file")
parser.add_argument("--model-bls",
                    dest="model_bls",
                    help="Path to model with branch length score")
parser.add_argument("--model-nobls",
                    dest="model_nobls",
                    help="Path to model without branch length score")
parser.add_argument("--output",
                    dest="output",
                    help="Output file name")

try:
    options = parser.parse_args()
except Exception, e:
    parser.print_help()

# redefine a functions for writing to stdout and stderr to save some writting
syserr = sys.stderr.write
sysout = sys.stdout.write


def main():
    """Main logic of the script"""
    #
    # Define mapping from column names to pretty names
    #
    columns_map = {'ContraScoreTargetSite': 'Accessibility',
                   'ID': "ID",
                   'MIRZABranchLengthScoreFill': "Conservation",
                   'MIRZAscore': 'MIRZAscore',
                   'distToBoundary': 'Distance to boundary',
                   'flanksG': 'Flanks G',
                   'flanksU': 'Flanks U',
                   'hybrid': "Hybrid",
                   'miRNA': 'miRNA',
                   'precise_type': 'Precise type',
                   'probability_with_bls': "Probability with conservation",
                   'probability_without_bls': "Probability without conservation",
                   'seed_beg': 'Seed start',
                   'seed_end': 'Seed end',
                   'type': 'Type'}
    #
    # Define columns order
    #
    columns_order = ['ID',
                     'miRNA',
                     'seed_beg',
                     'seed_end',
                     'type',
                     'precise_type',
                     'hybrid',
                     'flanksU',
                     'flanksG',
                     'ContraScoreTargetSite',
                     'distToBoundary',
                     'MIRZAscore',
                     'MIRZABranchLengthScoreFill',
                     'probability_without_bls',
                     'probability_with_bls']

    # read the models from pickle
    if options.verbose:
        syserr("Reading data file\n")
    df = pd.read_table(options.input)

    if options.verbose:
        syserr("Adding constant to data\n")
    df['const'] = 1.0


    model_bls =   sm.load(options.model_bls)
    model_nobls = sm.load(options.model_nobls)

    #
    # extract columns...
    #
    columns_bls = model_bls.params.keys().tolist()[1:]
    columns_nobls = model_nobls.params.keys().tolist()[1:]

    #
    # ...and predict and scale probabilities
    #
    if options.verbose:
        syserr("Adding probability to data\n - with BLS\n")
    df['probability_with_bls'] = scaled_logit_inverse(np.dot(df[['const'] + columns_bls].values, model_bls.params.values))
    if options.verbose:
        syserr(" - without BLS\n")
    df['probability_without_bls'] = scaled_logit_inverse(np.dot(df[['const'] + columns_nobls].values, model_nobls.params.values))

    #
    # reorder columns
    #
    if options.verbose:
        syserr("Reordering columns\n")
    df = df[columns_order]
    #
    # and rename columns
    #
    df.columns = [columns_map[col] for col in df.columns]

    if options.verbose:
        syserr("Saving file\n")
    df.to_csv(options.output, sep='\t', index=None, na_rep="NaN")
    if options.verbose:
        syserr("Done\n")

def mydot(x, y):
    """An implementation of the dot product
    hopefully without the seg fault error :)

    Args:
        x (@todo): @todo
        y (@todo): @todo

    Returns: @todo

    """
    return np.asarray([sum(starmap(mul,izip(i,y))) for i in x])

def scaled_logit_inverse(probability, scaling_factor=0.24):
    """This is inversed logit function with implemented
    scaling in the same time. It is equivalent of predicting
    probabilities with statsmodels and taking scale_prob after

    Args:
        probability (@todo): @todo

    Kwargs:
        scaling_factor (@todo): @todo

    Returns: @todo

    """
    return (scaling_factor*np.exp(probability))/(scaling_factor*np.exp(probability)+1.0)


def scale_prob(probability, scaling_factor=0.24):
    """ Functions that scales probability by:

        P' = P*c /(P*c + (1 - P))

    Args:
        probability (float or numpy array): probability for scaling

    Kwargs:
        scaling_factor (float): scaling factor to use

    Returns: float or numpy array

    """
    return (probability*scaling_factor) / (probability*scaling_factor + (1 - probability))


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
