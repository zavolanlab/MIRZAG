#!/usr/bin/env cwl-runner

cwlVersion: v1.0

class: CommandLineTool
baseCommand: rg_calculate_flanks_composition.py

hints:
  DockerRequirement:
    dockerPull: localhost:5000/zavolanlab/mirzag-scripts:1

doc: |

  Calculate flanks composition using provided coordinates.

arguments: ['--out', $(inputs.coords.basename).flanks]

inputs:

  seq:
    type: File
    doc: fasta with mRNA sequences
    inputBinding:
      prefix: --seq

  coords:
    type: File
    doc: file with target sites positions, miRNA and target gene ID
    inputBinding:
      prefix: --coords

  contextLen:
    type: ["null", int]
    default: 50
    doc: length of the context sequence serounding binding site
    inputBinding:
      prefix: --contextLen

  verbose:
    type: ["null", boolean]
    default: False
    doc: Be loud!
    inputBinding:
      prefix: --verbose


outputs:
  output:
    type: File
    outputBinding:
      glob: $(inputs.coords.basename).flanks
