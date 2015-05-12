#This script generates the jobs of the pipeline that will do the analysis on the splitted files.
import glob
import os
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
parser.add_argument("--group-id",
                    dest="group_id",
                    required=True,
                    help="Group Id")
parser.add_argument("--input-dir",
                    dest="input_dir",
                    required=True,
                    help="Input and output directory")

try:
    options = parser.parse_args()
except Exception, e:
    parser.print_help()

settings = ConfigObj(options.config).dict()

thisDir = os.getcwd()
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

for f in glob.glob(options.input_dir + "/*.fa"):
    input_name = os.path.splitext(f)[0]

    #
    # Calculate seed matches
    #
    seed_count_tuple = (os.path.join(pip_dir, "scripts/rg-count-miRNA-seeds.py"),
                  f,
                  settings['general']['seqs'],
                  input_name + ".seedcount")
    seed_count_command = "python %s --motifs %s --seqs %s --output %s --how TargetScan --coords --context 50" \
                          % seed_count_tuple
    seed_count_id = jobber.job(seed_count_command, {
        'name': 'SeedCount',
        'uniqueId': True})

    #
    # Filter duplicates if necessary
    #
    filter_duplicates_tuple = (os.path.join(pip_dir, "scripts/rg-filter-duplicates.py"),
                               input_name + ".seedcount",
                               os.path.join(options.input_dir, input_name + ".filterduplicates"),
                               settings['tasks']['FilterDuplicates']["split_by"],
                               settings['tasks']['FilterDuplicates']["index_after_split"]
                               )
    filter_duplicates_command = "python %s --coords %s --output %s --split-by %s --index-after-split %s -v" \
                                % filter_duplicates_tuple
    filter_duplicates_id = jobber.job(filter_duplicates_command, {
                                      'name': 'FilterDuplicates',
                                      'dependencies': [seed_count_id],
                                      'uniqueId': True})
    #
    # Create group for each file in order to calculate features
    #
    features_group_id = jobber.startGroup({'name': "Features_Group"})
    #
    # MIRZA
    #
    mirza_settings = settings['tasks']['CalculateMIRZA']
    mirza_script = 'scripts/rg-calculate-MIRZA.py'
    calculate_mirza_command = """python {script} \\
                    --out {output} \\
                    --seq {seqs} \\
                    --coords {input} \\
                    --motifs {motifs} \\
                    --contextLen {context} \\
                    --reforg {reforg} \\
                    --tree {tree} \\
                    --mln-dir {mlndir} \\
                    --threshold {threshold} \\
                    --onlyMIRZA {onlymirza} \\
                    --mirzabin {mirzabin} \\
                    -v
              """.format(**{'script': os.path.join(pip_dir, mirza_script),
                            'output': input_name + ".mirza",
                            'input': input_name + ".filterduplicates",
                            'seqs': settings['general']['seqs'],
                            'motifs': settings['general']['motifs'],
                            'context': mirza_settings.get('context_length', 50),
                            'reforg': mirza_settings.get('reference_organism', 'any/path'),
                            'tree': mirza_settings.get('phylogenetic_tree', 'any/path'),
                            'mlndir': mirza_settings.get('alignment_directory', 'any/path'),
                            'threshold': mirza_settings.get('threshold', 50),
                            'onlymirza': mirza_settings.get('run_only_MIRZA', "yes"),
                            'mirzabin': settings['general']['mirza_binary'],
                            })

    calculate_mirza_id = jobber.job(calculate_mirza_command, {
                                      'name': 'CalculateMIRZA',
                                      'dependencies': [filter_duplicates_id],
                                      'uniqueId': True})
    #
    # Contrafold
    #
    contrafold_settings = settings['tasks']['CalculateCONTRAfold']
    contrafold_script = 'scripts/rg-calculate-contrafold.py'
    contrafold_command = """python {script} \\
                                --out {output} \\
                                --seq {seqs} \\
                                --coords {input} \\
                                --contextLen_L {contextlen_l} \\
                                --contextLen_U {contextlen_u} \\
                                --context {context} \\
                                --contrabin {contrabin} \\
                                -v
                          """.format(**{'script': os.path.join(pip_dir, contrafold_script),
                                'output': input_name + ".contrafold",
                                'input': input_name + ".filterduplicates",
                                'seqs': settings['general']['seqs'],
                                'context': contrafold_settings.get('context', 50),
                                'contextlen_l': contrafold_settings.get('contextLen_L', 14),
                                'contextlen_u': contrafold_settings.get('contextLen_U', 0),
                                'contrabin': settings['general']['contrafold_binary']
                                 })
    calculate_contrafold_id = jobber.job(contrafold_command, {
                                      'name': 'CalculateCONTRAfold',
                                      'dependencies': [filter_duplicates_id],
                                      'uniqueId': True})
    #
    # Flanks
    #
    calculate_flanks_settings = settings['tasks']['CalculateFlanks']
    flanks_script = 'scripts/rg-calculate-flanks-composition.py'
    flanks_command = """python {script} \\
                            --out {output} \\
                            --seq {seqs} \\
                            --coords {input} \\
                            --contextLen {context} \\
                            -v
                      """.format(**{'script': os.path.join(pip_dir, flanks_script),
                                   'output': input_name + ".flanks",
                                   'input': input_name + ".filterduplicates",
                                   'seqs': settings['general']['seqs'],
                                   'context': calculate_flanks_settings.get('context_length', 50),
                                 })
    calculate_flanks_id = jobber.job(flanks_command, {
                                      'name': 'CalculateFlanks',
                                      'dependencies': [filter_duplicates_id],
                                      'uniqueId': True})

    #
    # Distance
    #
    calculate_distance_settings = settings['tasks']['CalculateDistance']
    distance_script = 'scripts/rg-calculate-distance.py'
    distance_command = """python {script} \\
                            --out {output} \\
                            --seq {seqs} \\
                            --coords {input} \\
                            -v
                      """.format(**{'script': os.path.join(pip_dir, distance_script),
                                     'output': input_name + ".distance",
                                     'input': input_name + ".filterduplicates",
                                     'seqs': settings['general']['seqs'],
                                    })
    calculate_distance_id = jobber.job(distance_command, {
                                      'name': 'CalculateDistance',
                                      'dependencies': [filter_duplicates_id],
                                      'uniqueId': True})
    jobber.endGroup()
    #
    # Merge and add probabilities
    #
    merge_script = 'scripts/rg-merge-results-and-add-probability.py'
    merge_settings = settings['tasks']['MergeAndAddProbability']
    merge_inputs = ",".join([input_name + ".contrafold",
                             input_name + ".mirza",
                             input_name + ".flanks",
                             input_name + ".distance"])
    merge_command = """python {script} \\
                    --output {output} \\
                    --inputs {inputs} \\
                    --coords {coords} \\
                    --model-bls {model_bls} \\
                    --model-nobls {model_nobls} \\
                    --only-mirza {onlymirza} \\
                    -v
              """.format(**{'script': os.path.join(pip_dir, merge_script),
                                  'output': input_name + ".probabilities",
                                  'inputs': merge_inputs,
                                  'coords': input_name + ".filterduplicates",
                                  'model_bls':   settings['general']['model_with_bls'],
                                  'model_nobls': settings['general']['model_without_bls'],
                                  'onlymirza': merge_settings.get('run_only_MIRZA', 'yes')
                                 })
    merge_id = jobber.job(merge_command, {
                                      'name': 'MergeAndAddProbability',
                                      'dependencies': [features_group_id],
                                      'uniqueId': True})

    #
    # Calculate score
    #
    calculate_score_script = 'scripts/rg-calculate-per-gene-score.py'
    calculate_score_settings = settings['tasks']['CalculateScore']
    calculate_score_command = """python {script} \\
                    --output {output} \\
                    --input {input} \\
                    --threshold {threshold} \\
                    --split-by "{split_by}" \\
                    --colum {column} \\
                    -v
              """.format(**{'script': os.path.join(pip_dir, calculate_score_script),
                                  'output': input_name + ".score",
                                  'input': input_name + ".probabilities",
                                  'threshold': calculate_score_settings.get('threshold', 0.12),
                                  'split_by':  calculate_score_settings.get('split_by', "NOTHING"),
                                  'column':    calculate_score_settings.get('index_after_split', 0),
                                  })
    merge_id = jobber.job(calculate_score_command, {
                                      'name': 'CalculateScore',
                                      'dependencies': [merge_id],
                                      'uniqueId': True})

jobber.endGroup()
jobber.launch(options.group_id)
