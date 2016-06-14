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
import csv
import time
import gzip
import pandas as pd
import numpy as np
import statsmodels.api as sm
from argparse import ArgumentParser

parser = ArgumentParser(description=__doc__)
parser.add_argument("-v",
                    "--verbose",
                    dest="verbose",
                    action="store_true",
                    default=False,
                    help="Be loud!")
parser.add_argument("--inputs",
                    dest="inputs",
                    help="Coma-separated list of input paths")
parser.add_argument("--coords",
                    dest="coords",
                    help="Cooridinate file  used in the beginning")
parser.add_argument("--model-bls",
                    dest="model_bls",
                    help="Path to model with branch length score")
parser.add_argument("--model-nobls",
                    dest="model_nobls",
                    help="Path to model without branch length score")
parser.add_argument("--output",
                    dest="output",
                    help="Output file name")
parser.add_argument("--only-mirza",
                    choices=('yes', 'no'),
                    dest="only_mirza",
                    help="Calculate only MIRZA and DON'T calculate MIRZA BLS")


# redefine a functions for writing to stdout and stderr to save some writting
syserr = sys.stderr.write
sysout = sys.stdout.write

class EmptyDataException(Exception): pass


def main(options):
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
                     'flanksU',
                     'flanksG',
                     'ContraScoreTargetSite',
                     'distToBoundary',
                     'MIRZAscore',
                     'MIRZABranchLengthScoreFill',
                     'probability_without_bls',
                     'probability_with_bls']
    #
    # read the data and clean it a little
    #
    accessibility, mirza, flanks, distance = options.inputs.split(",")

    if options.verbose:
        syserr("Preparing coordinates\n")
    with gzip.open(options.coords) as gcor:
        data = {",".join(row[:-1]): {} for row in csv.reader(gcor, delimiter='\t')}

    if options.verbose:
        syserr(" - Adding Accessibility\n")
    for (id, ascore) in csv.reader(gzip.open(accessibility), delimiter='\t'):
        # print id, ascore
        if id in data:
            data[id].update({'ContraScoreTargetSite': ascore})

    if options.verbose:
        syserr(" - Adding MIRZA\n")
    for (id, score, bls) in csv.reader(gzip.open(mirza), delimiter='\t'):
        if id in data:
            data[id].update({'MIRZAscore': score,
                             'MIRZABranchLengthScoreFill': bls})

    if options.verbose:
        syserr(" - Adding Distance\n")
    for (id, dist) in csv.reader(gzip.open(distance), delimiter='\t'):
        if id in data:
            data[id].update({'distToBoundary': dist})

    if options.verbose:
        syserr(" - Adding Flanks\n")
    for (id, flanksG, flanksA, flanksC, flanksU) in csv.reader(gzip.open(flanks), delimiter='\t'):
        if id in data:
            data[id].update({'flanksG': flanksG,
                             'flanksU': flanksU})
    #
    # if options.verbose:
    #     syserr("Converting dictionary to data\n")
    # # for the testing reasons lets save the dataframe
    # if options.verbose:
    #     syserr(" - Saving data\n")
    # save(data, os.path.join(options.output_dir, "chunk_%i.bin" % (counter + 1)))
    data = pd.DataFrame(data).T.replace("NA", np.nan)
    gc.collect()

    if options.verbose:
        syserr("Adding additional columns:\n")
    if options.verbose:
        syserr(" - id\n")
    data['ID'] = [i.split(',')[0] for i in data.index]
    if options.verbose:
        syserr(" - miRNA\n")
    data['miRNA'] = [i.split(',')[1] for i in data.index]
    if options.verbose:
        syserr(" - seed beginning\n")
    data['seed_beg'] = [i.split(',')[2] for i in data.index]
    if options.verbose:
        syserr("seed end\n")
    data['seed_end'] = [i.split(',')[3] for i in data.index]

    if options.verbose:
        syserr("Transforming MIRZA to log space\n")

    if len(data) == 0:
        with gzip.open(options.output, 'wb') as handler:
            pass
        raise EmptyDataException("Empty DataFrame - exiting")

    data['MIRZAscore'] = np.log(data['MIRZAscore'].astype(np.float))


    if options.verbose:
        syserr("Adding constant to data\n")
    data['const'] = 1.0
    # import ipdb; ipdb.set_trace()

    model_bls = pd.read_pickle(options.model_bls)
    model_nobls = pd.read_pickle(options.model_nobls)

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
    if options.only_mirza == 'no':
        if options.verbose:
            syserr("    > doing dot product\n")
        mydot = dot_product(data[['const'] + columns_bls].astype(np.float).values, model_bls.params.values)
        if options.verbose:
            syserr("    > scaling probability\n")
        data['probability_with_bls'] = scaled_logit_inverse(mydot)
    else:
        if options.verbose:
            syserr("    > adding only NaN because only-mirza = yes\n")
        data['probability_with_bls'] = np.nan

    if options.verbose:
        syserr(" - without BLS\n")
    if options.verbose:
        syserr("    > doing dot product\n")
    mydot = dot_product(data[['const'] + columns_nobls].astype(np.float).values, model_nobls.params.values)
    if options.verbose:
        syserr("    > scaling probability\n")
    data['probability_without_bls'] = scaled_logit_inverse(mydot)

    #
    # reorder columns
    #
    if options.verbose:
        syserr("Reordering columns\n")
    data = data[columns_order]
    #
    # and rename columns
    #
    data.columns = [columns_map[col] for col in data.columns]

    if options.verbose:
        syserr("Saving file\n")
    import ipdb; ipdb.set_trace()
    with gzip.open(options.output, 'wb') as handler:
        data.to_csv(handler, sep='\t', index=None, na_rep="NaN", header=None)


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


def dot_product(x, y):
    """A dot product of the matrix and the vecotr

    Args:
        x (np.array): numpy array of vectors of length y
        y (np.array): one dimensional numpy array

    Returns: np.array

    """
    return np.asarray([np.dot(i, y) for i in x])

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
        try:
            main(options)
        except EmptyDataException, e:
            syserr(str(e) + "\n")
        if options.verbose:
            syserr("### Successfully finished in %i seconds, on %s ###\n" % (time.time() - start_time, time.strftime("%d-%m-%Y at %H:%M:%S")))
    except KeyboardInterrupt:
        syserr("Interrupted by user after %i seconds!\n" % (time.time() - start_time))
        sys.exit(-1)
