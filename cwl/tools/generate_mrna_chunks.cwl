#!/usr/bin/env cwl-runner

cwlVersion: v1.0

class: CommandLineTool
baseCommand: rg_generate_utr_chunks.py

doc: |

  Take mRNAs and generate fragments by sliding window that can be
  fed into MIRZA


inputs:

  input:
    type: File
    doc: Input file in FASTA format containing mRNA sequences
    inputBinding:
      prefix: --input

  part_size:
    type: ["null", int]
    default: 40000
    doc: Number of sequences per part, defaults to 40000
    inputBinding:
      prefix: --part-size

  window_size:
    type: ["null", int]
    default: 50
    doc: Length of the window for MIRZA, defaults to 50
    inputBinding:
      prefix: --window-size

  slide_size:
    type: ["null", int]
    default: 20
    doc: Size of the window slide, defaults to 20
    inputBinding:
      prefix: --slide-size

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
      glob: "*.mrna.fa"
