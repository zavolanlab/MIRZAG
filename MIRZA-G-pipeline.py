import os
import random
import sys

from configobj import ConfigObj
from argparse import ArgumentParser, RawTextHelpFormatter

parser = ArgumentParser(description=__doc__, formatter_class=RawTextHelpFormatter)
parser.add_argument("-v",
                    "--verbose",
                    dest="verbose",
                    action="store_true",
                    default=False,
                    help="Be loud!")
parser.add_argument("--config",
                    dest="config",
                    required=True,
                    help="Config file")

try:
    options = parser.parse_args()
except Exception, e:
    parser.print_help()

settings = ConfigObj(options.config).dict()

working_directory = os.getcwd()
pipeline_directory = os.path.dirname(os.path.abspath(__file__))
output_directory = os.path.join(working_directory, "output")
os.makedirs(output_directory)
jobber_path = settings['general']['jobber_path']
sys.path.append(jobber_path)
from jobber import JobClient

jobber = JobClient.Jobber()

#Create a group for whole pipeline. The module "Python" will be inherited by all jobs that are in this group,
# so we don't need to define it for each job that calls a python script
pipeline_id = jobber.startGroup({'name': "MIRZA-G",
                                 'options': [['module', "Python"],
                                             ['module', "GCC"],
                                             ['module', "CONTRAfold"]],
                                 'executer': settings['general'].get('executer', 'drmaa')})

#First step is to split the file
split_command = "python %s --input %s --output-dir %s" % (os.path.join(pipeline_directory, "scripts/rg-prepare-mirnas-for-mirza-and-split.py"),
                                                          settings['general']['motifs'],
                                                          output_directory)
split_files_id = jobber.job(split_command, {'name': "SplitMiRNAs"})


#We create a group where the jobs to analyse the splitted files will be put into
analyse_files_id = jobber.startGroup({'name': "analysis",
                                    'dependencies': [split_files_id]})

#We call the script that will generate the jobs that will analyse the split files. We pass the id of the group
#and the folder where the script will find the splitted files.
analysis_tuple = (os.path.join(pipeline_directory, "run-analysis.py"),
                  output_directory,
                  analyse_files_id,
                  os.path.abspath(options.config))
analysis_command = "python %s --input-dir %s --group-id %s --config %s -v" % analysis_tuple
jobber.job(analysis_command, {'name': "createJobs"})


jobber.endGroup()

# We merge the files into our result file after analysis finishes
final_merge_command = "zcat {output_dir}/*.score > {cwd}/mirza_g_results.tab".format(output_dir=output_directory,
                                                                                     cwd=working_directory)
jobber.job(final_merge_command, {'name': "MergeScore",
                                 'dependencies': [analyse_files_id]})


jobber.endGroup()

# Before launching we print the command to stop the pipeline
print "In order to stop the pipeline run a command:"
print "python %s/jobber_server.py -command delete -jobId %i" % (jobber_path, pipeline_id)

#You need to always launch, otherwise jobs wont get executed.
jobber.launch(pipeline_id)
