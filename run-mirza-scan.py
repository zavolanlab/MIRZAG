#This script generates the jobs of the pipeline that will do the analysis on the splitted files.
import glob
import os
import sys
from configobj import ConfigObj
from jinja2 import Template
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
parser.add_argument("--group-id",
                    dest="group_id",
                    required=True,
                    help="Group Id")
parser.add_argument("--input-dir",
                    dest="input_dir",
                    required=True,
                    help="Input and output directory")
parser.add_argument("--working-dir",
                    dest="working_dir",
                    required=True,
                    help="Working directory of the pipeline. Required because this file is launched from ~")

try:
    options = parser.parse_args()
except Exception, e:
    parser.print_help()

settings = ConfigObj(options.config).dict()

thisDir = options.working_dir
pip_dir = os.path.dirname(os.path.abspath(__file__))
jobberDir = settings['general']['jobber_path']
sys.path.append(jobberDir)
from jobber import JobClient

jobber = JobClient.Jobber()

#We can call "extendGroup" if we want to create jobs into an already existing group. Don't forget to call "endGroup" and "launch"
#after you're done
jobber.extendGroup(options.group_id)

#We assume that each file to be analyzed ends with .seqs. Its important to always distinguish input files from any intermediary
#files in case we need to restart the jobs. We should make the jobs unique to prevent duplication of jobs in case this script
#is run multiple times

# template = settings['general']['template']
template = os.path.join(pip_dir, "scripts/template.sh")
with open(template) as tmpl:
    template = Template(tmpl.read())

scan_group = jobber.startGroup({'name': 'CalculateCoordinates'})
for f in glob.glob(options.input_dir + "/*.fa"):
    input_name = os.path.splitext(f)[0]

    #
    # Calculate seed matches
    #
    scan_settings = settings['tasks']['ScanWithMIRZA']
    scan_script = 'scripts/rg-extract-data-from-mirza-output.py'
    scan_command = """{mirza_bin} {expressions} {mrnas} {mirnas} 50 noupdate | python {script} \\
                                --seqs {seqs} \\
                                --output {output} \\
                                --context {context} \\
                                --threshold {threshold} \\
                                -v
                  """
    #
    # If there is template use it for command
    #
    if settings['general'].get('executer', 'drmaa') == 'drmaa':
        #
        # Copy files by default to the tmp directory
        #
        copy_dir = "$TMPDIR"
        copy_files = {f: 'input.fa',
                      settings['general']['seqs']: 'seqs.fa',
                      os.path.join(options.input_dir, "mirnas.expression"): 'expressions',
                      settings['general']['motifs']: 'motifs.fa'}
        moveback = {'output': input_name + ".mirzascan"}

        seed_command_rendered = template.render(modules=scan_settings.get('modules', None),
                                           command=scan_command,
                                           copy=copy_files,
                                           moveback=moveback,
                                           copydir=copy_dir)
        scan_command = str(seed_command_rendered).format(**{
                                       'mirza_bin': settings['general']['mirza_binary'],
                                       'mrnas': 'input.fa',
                                       'mirnas': 'motifs.fa',
                                       'expressions': 'expressions',
                                       'script': os.path.join(pip_dir, scan_script),
                                       'seqs': 'seqs.fa',
                                       'output': 'output',
                                       'threshold': scan_settings.get('threshold', 50),
                                       'context': scan_settings.get('context', 50)})
    else:
        scan_command = str(seed_command_rendered).format(**{
                                       'mirza_bin': settings['general']['mirza_binary'],
                                       'mrnas': f,
                                       'mirnas': settings['general']['motifs'],
                                       'expressions': os.path.join(options.input_dir, "mirnas.expression"),
                                       'script': os.path.join(pip_dir, scan_script),
                                       'seqs': settings['general']['seqs'],
                                       'output': input_name + ".mirzascan",
                                       'threshold': scan_settings.get('threshold', 50),
                                       'context': scan_settings.get('context', 50)})
    scan_id = jobber.job(scan_command,
                               {'name': 'ScanWithMIRZA',
                                'uniqueId': True,
                                'options': [('q', scan_settings.get('queue', 'short.q')),
                                            ('l', "h_vmem=%s" % scan_settings.get('mem_req', '2G'))]
                                })
jobber.endGroup()


# We merge the files into our result file after analysis finishes
merge_command = "cat {output_dir}/*.mirzascan > {output_dir}/scan_results.tab".format(output_dir=options.input_dir)
merge_id = jobber.job(merge_command, {'name': "MergeScan",
                                   'dependencies': [scan_group]})

# Filter results

filter_results_settings = settings['tasks']['FilterScan']
filter_results_script = 'scripts/rg-filter-duplicates-from-scan.py'
filter_results_command = """python {script} \\
                                --coords {input} \\
                                --output {output} \\
                                --split-by \"{split_by}\" \\
                                --index-after-split {index_after_split} \\
                                -v
              """
#
# If there is template use it for command
#
if settings['general'].get('executer', 'drmaa') == 'drmaa':
    #
    # Copy files by default to the tmp directory
    #
    copy_dir = "$TMPDIR"
    copy_files = {os.path.join(options.input_dir, "scan_results.tab"): 'input'}
    moveback = {'output': os.path.join(options.input_dir, "scan_result.filtered")}

    filter_results_command_rendered = template.render(modules=filter_results_settings.get('modules', None),
                                       command=filter_results_command,
                                       copy=copy_files,
                                       moveback=moveback,
                                       copydir=copy_dir)
    filter_results_command = str(filter_results_command_rendered).format(**{'script': os.path.join(pip_dir, filter_results_script),
                                     'input': 'input',
                                     'output': 'output',
                                     'index_after_split': settings['general'].get('index_after_split'),
                                     'split_by': settings['general'].get('split_by', "NONE"),
                                   })
else:
    filter_results_command = str(filter_results_command).format(**{'script': os.path.join(pip_dir, filter_results_script),
                                     'input': os.path.join(options.input_dir, "scan_results.tab"),
                                     'output': os.path.join(options.input_dir, "scan_result.filtered"),
                                     'index_after_split': settings['general'].get('index_after_split'),
                                     'split_by': settings['general'].get('split_by', "NONE"),
                                   })
filter_results_id = jobber.job(filter_results_command,
                           {'name': 'FilterScan',
                            'uniqueId': True,
                             'dependencies': [merge_id],
                            'options': [('q', filter_results_settings.get('queue', 'short.q')),
                                        ('l', "h_vmem=%s" % filter_results_settings.get('mem_req', '2G'))]
                            })

jobber.endGroup()
jobber.launch(options.group_id)
