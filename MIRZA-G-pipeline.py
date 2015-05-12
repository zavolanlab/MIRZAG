#This example shows you how to make a pipeline where we split a file
#and create jobs based on the number of files we created. After the jobs finish,
#we merge the results back. Note that in this example, we do not guarantee the order
#of the results. So the merging can be arbitrary. If you need order, you would need to ensure this yourself,
#e.g. by sorting the files by filename before reading them.
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

thisDir = os.getcwd()
pip_dir = os.path.dirname(os.path.abspath(__file__))
output_dir = os.path.join(thisDir, "output")
os.makedirs(output_dir)
jobberDir = settings['general']['jobber_path']
sys.path.append(jobberDir)
from jobber import JobClient

jobber = JobClient.Jobber()

#Create a group for whole pipeline. The module "Python" will be inherited by all jobs that are in this group,
# so we don't need to define it for each job that calls a python script
pipelineId = jobber.startGroup({
    'name': "Pipeline",
    'options': [
        ['module', "Python"],
        ['module', "GCC"],
        ['module', "CONTRAfold"],
    ]
})

#First step is to split the file
split_command = "python %s --input %s --output-dir %s" % (os.path.join(pip_dir, "scripts/rg-prepare-mirnas-for-mirza-and-split.py"),
                                                          settings['general']['motifs'],
                                                          output_dir)
splitFilesId = jobber.job(split_command, {
    'name': "SplitMiRNAs"
})

#We create a group where the jobs to analyse the splitted files will be put into

analyseFilesId = jobber.startGroup({
    'name': "analysis",
    'dependencies': [splitFilesId]
})

#We call the script that will generate the jobs that will analyse the splitted files. We pass the id of the group
#and the folder where the script will find the splitted files.
analysis_tuple = (os.path.join(pip_dir, "run-analysis.py"),
                  output_dir,
                  analyseFilesId,
                  os.path.abspath(options.config))
analysis_command = "python %s --input-dir %s --group-id %s --config %s -v" % analysis_tuple
jobber.job(analysis_command, {
    'name': "createJobs"
})


jobber.endGroup()

# We merge the files into our result file after analysis finishes
jobber.job("zcat " + output_dir + "/*.score > " + thisDir + "/mirza_results.tab", {
    'name': "MergeScore",
    'dependencies': [analyseFilesId]
})


jobber.endGroup()
#You need to always launch, otherwise jobs wont get executed.
jobber.launch(pipelineId)
