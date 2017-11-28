#!/usr/bin/env cwl-runner

cwlVersion: v1.0

class: CommandLineTool
baseCommand: rg_extract_data_from_mirza_output.py

doc: |

  Take MIRZA output and arrange it in a proper way

arguments: ['--output', $(inputs.mirza.basename).tsv]

inputs:

  mirza:
    type: File
    doc: Input file from MIRZA
    inputBinding:
      prefix: --input

  mrna:
    type: File
    doc: mRNA sequences in FASTA format
    inputBinding:
      prefix: --seqs

  threshold:
    type: ["null", float]
    default: 50
    doc: Threshold for the score
    inputBinding:
      prefix: --threshold

  context:
    type: ["null", int]
    default: 50
    doc: Context for sequence to print
    inputBinding:
      prefix: --context

  verbose:
    type: ["null", boolean]
    default: False
    doc: Be loud!
    inputBinding:
      prefix: --verbose

stdout: $(inputs.mirza.basename).tsv

outputs:
  output:
    type: File
    outputBinding:
      glob: $(inputs.mirza.basename).tsv
