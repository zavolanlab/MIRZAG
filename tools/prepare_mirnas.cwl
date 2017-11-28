#!/usr/bin/env cwl-runner

cwlVersion: v1.0

class: CommandLineTool
baseCommand: rg_prepare_mirnas_for_mirza_and_split.py

doc: |

  Prepare miRNA fasta file for MIRZA i.e. for each miRNA sequence
  check if it is 21 nucleotide long and if not eliminate it. It
  also replaces all u or U into T. In the same time it splits miRNAs
  into separate files.


inputs:

  input:
    type: File
    doc: Input miRNA file in FASTA format.
    inputBinding:
      prefix: --input

  verbose:
    type: ["null", boolean]
    default: False
    doc: Be loud!
    inputBinding:
      prefix: --verbose


outputs:
  output:
    type: File[]
    outputBinding:
      glob: "*.mirna.fa"
