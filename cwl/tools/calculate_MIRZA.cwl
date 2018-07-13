#!/usr/bin/env cwl-runner

cwlVersion: v1.0

class: CommandLineTool
baseCommand: rg_calculate_MIRZA.py

hints:
  DockerRequirement:
    dockerPull: zavolab/mirzag-scripts:1

doc: |

  Calculate MIRZA interaction energy and MIRZA-based Branch Length Score
  (conservation) for the provided coordinates of putative targets.

arguments: ['-o', $(inputs.coords.basename).mirza]

inputs:

  seq:
    type: File
    doc: Fasta with mRNA sequences
    inputBinding:
      prefix: --seq

  motifs:
    type: File
    doc: Fasta file with miRNA sequences
    inputBinding:
      prefix: --motifs

  coords:
    type: File
    doc: File with target sites positions, miRNA and target gene ID
    inputBinding:
      prefix: --coords

  tree:
    type: File
    doc: Phylogenetic tree of the species used in alignment file
    inputBinding:
      prefix: --tree

  msa:
    type: ['null', File]
    default: null
    doc: Archive (tar.gz) with multiple sequence alignment files
    inputBinding:
      prefix: --msa

  onlymirza:
    type:
    - "null"
    - type: enum
      symbols: ['yes', 'no']
    default: no
    doc: Calculate only MIRZA score for given coordinates
    inputBinding:
      prefix: --onlyMIRZA

  threshold:
    type: ["null", float]
    default: 20.0
    doc: Threshold for MIRZA score
    inputBinding:
      prefix: --threshold

  contextLen:
    type: ["null", int]
    default: 50
    doc: Length of the context sequence serounding binding site
    inputBinding:
      prefix: --contextLen

  mirzabin:
    type: ["null", string]
    default: MIRZA
    doc: Path to the MIRZA binary
    inputBinding:
      prefix: --mirzabin

  reforg:
    type: ["null", string]
    default: hg19
    doc: Reference organism to which alignments are performed (default=hg19)
    inputBinding:
      prefix: --reforg

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
      glob: $(inputs.coords.basename).mirza
