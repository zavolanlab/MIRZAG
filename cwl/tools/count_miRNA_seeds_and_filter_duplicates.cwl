#!/usr/bin/env cwl-runner

cwlVersion: v1.0

class: CommandLineTool
baseCommand: rg_count_miRNA_seeds_and_filter_duplicates.py

hints:
  DockerRequirement:
    dockerPull: localhost:5000/zavolanlab/mirzag-scripts:1

doc: |

  Count miRNA seed in the provided sequences according to chosen definition

arguments: ['--output', $(inputs.motifs.basename).seedcount]

inputs:

  motifs:
    type: File
    doc: miRNA/siRNA sequences to use when scanning
    inputBinding:
      prefix: --motifs

  seqs:
    type: File
    doc: "Sequences for scaning eg. 3' UTRs"
    inputBinding:
      prefix: --seqs

  how:
    type:
    - "null"
    - type: enum
      symbols: ['ElMMo', 'TargetScan', '6-mer']
    default: TargetScan
    doc: What definition for seed to use
    inputBinding:
      prefix: --how

  context:
    type: ["null", int]
    default: 50
    doc: Context for sequence to print, defaults to 50
    inputBinding:
      prefix: --context

  split_by:
    type: string
    doc: Split id by the string
    inputBinding:
      prefix: --split-by

  index_after_split:
    type: int
    doc: After split take this column as new id, 0 based
    inputBinding:
      prefix: --index-after-split

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
      glob: $(inputs.motifs.basename).seedcount
