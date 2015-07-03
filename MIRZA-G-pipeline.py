import os
import sys
import shutil
import random
import traceback
from errno import EEXIST

from configobj import ConfigObj
from argparse import ArgumentParser, RawTextHelpFormatter


parser = ArgumentParser(description=__doc__, formatter_class=RawTextHelpFormatter)
subparsers = parser.add_subparsers(help="Commands", dest='command')
run_parser = subparsers.add_parser("run", help="Run a pipeline")
run_parser.add_argument("-v",
                    "--verbose",
                    dest="verbose",
                    action="store_true",
                    default=False,
                    help="Be loud!")
run_parser.add_argument("--config",
                    dest="config",
                    required=True,
                    help="Config file")
run_parser.add_argument("--name-suffix",
                    dest="name_suffix",
                    default="test_run",
                    help="Suffix to add to pipeline name in order to easily differentiate between different run, defaults to test_run")
run_parser.add_argument("--protocol",
                    dest="protocol",
                    default="seed",
                    choices=("seed", "scan"),
                    help="Protocol of MIRZA-G, defaults to seed")
run_parser.add_argument("--calulate-bls",
                    dest="calulate_bls",
                    action="store_true",
                    default=False,
                    help="Calculate Branch Length Score (conservation)")

clean_parser = subparsers.add_parser("clean", help="Clean after previous run")
clean_parser.add_argument("-v",
                    "--verbose",
                    dest="verbose",
                    action="store_true",
                    default=False,
                    help="Be loud!")
clean_parser.add_argument("-y",
                        "--yes",
                        dest="yes",
                        action="store_true",
                        default=False,
                        help="Force deletion of files.")

try:
    options = parser.parse_args()
except Exception, e:
    parser.print_help()


def mkdir_p(path_to_dir):
    """Make directory with subdirectories"""
    try:
        os.makedirs(path_to_dir)
    except OSError as e: # Python >2.5
        if e.errno == EEXIST and os.path.isdir(path_to_dir):
            message = "Output directory %s (and files) exist.\nYou need to clean before you proceed. Run:\n" + \
            "python MIRZA-G_pipeline clean \n"
            sys.stderr.write(message % path_to_dir)
            sys.exit()
        else:
            raise e

# redefine a functions for writing to stdout and stderr to save some writting
syserr = sys.stderr.write
sysout = sys.stdout.write


working_directory = os.getcwd()
pipeline_directory = os.path.dirname(os.path.abspath(__file__))
output_directory = os.path.join(working_directory, "output")


if options.command == 'clean':
    try:
        if options.yes:
            is_sure = "yes"
        else:
            is_sure = raw_input("Do you really want to delete previous run (yes/no)?:  ")
        if is_sure.upper().startswith("Y"):
            try:
                shutil.rmtree(output_directory)
            except OSError:
                if options.verbose:
                    syserr(" -> no such file or directory: %s\n" % output_directory)
            files_to_delete = ["mirza_g_results_scan.tab",
                               "mirza_g_results_seed.tab"]
            for f in files_to_delete:
                if options.verbose:
                    syserr("Removing %s\n" % os.path.join(working_directory, f))
                try:
                    os.remove(os.path.join(working_directory, f))
                except OSError, e:
                    if options.verbose:
                        syserr(" -> no such file or directory: %s\n" % os.path.join(working_directory, f))
            if options.verbose:
                syserr("All output files and directories were cleaned\n")
    except Exception as e:
        syserr(traceback.format_exc())
    finally:
        sys.exit()

settings = ConfigObj(options.config).dict()
mkdir_p(output_directory)
if options.protocol == "scan":
    mkdir_p(os.path.join(output_directory, "MIRZAscan"))
jobber_path = settings['general']['jobber_path']
sys.path.append(jobber_path)
from jobber import JobClient

jobber = JobClient.Jobber()

#Create a group for whole pipeline. The module "Python" will be inherited by all jobs that are in this group,
# so we don't need to define it for each job that calls a python script
pipeline_id = jobber.startGroup({'name': "MIRZA-G_%s" % options.name_suffix,
                                 'options': [['module', "Python"],
                                             ['module', "GCC"],
                                             ['module', "CONTRAfold"]],
                                 'executer': settings['general'].get('executer', 'drmaa')})

if options.protocol == "seed":

    #First step is to split the file
    split_command = "python %s --input %s --output-dir %s" % (os.path.join(pipeline_directory, "scripts/rg-prepare-mirnas-for-mirza-and-split.py"),
                                                              settings['general']['motifs'],
                                                              output_directory)
    split_files_id = jobber.job(split_command, {'name': "SplitMiRNAs"})

if options.protocol == "scan":
    #First step is to split the file
    make_chunks = "python %s --input %s --output-dir %s -v" % (os.path.join(pipeline_directory, "scripts/rg-generate-utr-chunks.py"),
                                                                 settings['general']['seqs'],
                                                                 os.path.join(output_directory, "MIRZAscan"))
    make_chunks_id = jobber.job(make_chunks, {'name': "GenerateChunks"})

    # generate expressions
    gen_expr_tup = (settings['general']['motifs'],
                    os.path.join(output_directory, "MIRZAscan/mirnas.expression"))
    gen_expr_command = """cat %s | ruby -ne 'puts "#{$_.rstrip()[1..-1]}\t1" if $_.start_with?(">")' > %s""" % gen_expr_tup
    gen_expressions_id = jobber.job(gen_expr_command, {'name': "GenerateExpressions"})

    # We create a group where the jobs to analyse the splitted files will be put into and MIRZA will be calculated
    mirza_scan_id = jobber.startGroup({'name': "MIRZAscan",
                                       'dependencies': [make_chunks_id, gen_expressions_id]})

    #We call the script that will generate the jobs that will analyse the split files. We pass the id of the group
    #and the folder where the script will find the splitted files.
    scan_tuple = (os.path.join(pipeline_directory, "run-mirza-scan.py"),
                      os.path.join(output_directory, "MIRZAscan"),
                      mirza_scan_id,
                      os.path.abspath(options.config),
                      working_directory)
    scan_command = "python %s --input-dir %s --group-id %s --config %s -v --working-dir %s" % scan_tuple
    jobber.job(scan_command, {'name': "createScanJobs"})


    jobber.endGroup()

    split_command = """awk '{print | "gzip > %s/"$2".seedcount"}' %s""" % (output_directory,
                                                                    os.path.join(output_directory, "MIRZAscan/scan_result.filtered"))

    split_files_id = jobber.job(split_command, {'name': "SplitCoords",
                                                'dependencies': [mirza_scan_id]})


#We create a group where the jobs to analyse the splitted files will be put into
analyse_files_id = jobber.startGroup({'name': "Analysis",
                                    'dependencies': [split_files_id]})

#We call the script that will generate the jobs that will analyse the split files. We pass the id of the group
#and the folder where the script will find the splitted files.
analysis_tuple = (os.path.join(pipeline_directory, "run-analysis.py"),
                  output_directory,
                  analyse_files_id,
                  os.path.abspath(options.config),
                  working_directory,
                  options.protocol)
analysis_command = "python %s --input-dir %s --group-id %s --config %s -v --working-dir %s --protocol %s" % analysis_tuple
jobber.job(analysis_command, {'name': "createJobs"})


jobber.endGroup()

# We merge the files into our result file after analysis finishes
final_merge_command = "zcat {output_dir}/*.score > {cwd}/mirza_g_results_{protocol}.tab".format(output_dir=output_directory,
                                                                                     cwd=working_directory,
                                                                                                protocol=options.protocol)
jobber.job(final_merge_command, {'name': "MergeResults",
                                 'dependencies': [analyse_files_id]})


jobber.endGroup()

# Before launching we print the command to stop the pipeline
print "In order to stop the pipeline run a command:"
print "python %s/jobber_server.py -command delete -jobId %i" % (jobber_path, pipeline_id)

#You need to always launch, otherwise jobs wont get executed.
jobber.launch(pipeline_id)
