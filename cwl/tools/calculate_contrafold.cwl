#!/usr/bin/env cwl-runner

cwlVersion: v1.0

class: CommandLineTool
baseCommand: rg_calculate_contrafold.py

hints:
  DockerRequirement:
    dockerPull: localhost:5000/zavolanlab/mirzag-scripts:1

doc: |

  The `CONTRAfold <http://contra.stanford.edu/contrafold/>`_ algorithm is used
  to calculate accessibility of target site for miRNA.

arguments: ['--out', $(inputs.coords.basename).contrafold]

inputs:

  seq:
    type: File
    doc: Fasta file with mRNA sequences
    inputBinding:
      prefix: --seq

  coords:
    type: File
    doc: file with target sites positions, miRNA and target gene ID
    inputBinding:
      prefix: --coords

  contextLen_L:
    type: ["null", int]
    default: 0
    doc: length of the context sequence downstream binding site to be unwinded
    inputBinding:
      prefix: --contextLen_L

  contextLen_U:
    type: ["null", int]
    default: 0
    doc: length of the context sequence upstream binding site to be unwinded
    inputBinding:
      prefix: --contextLen_U

  context:
    type: ["null", int]
    default: 50
    doc: length of the context of the seed to be checked
    inputBinding:
      prefix: --context

  contrabin:
    type: ["null", string]
    default: contrafold
    doc: Path to CONTRAfold binary
    inputBinding:
      prefix: --contrabin

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
      glob: $(inputs.coords.basename).contrafold
