#!/usr/bin/env cwl-runner

cwlVersion: v1.0

class: CommandLineTool
baseCommand: rg_merge_results_add_probability_and_calculate_per_gene_score.py

doc: |

  Merge all results into one features table

arguments: ['--output', $(inputs.coords.basename).score]

inputs:

  coords:
    type: File
    doc: File with target sites positions, miRNA and target gene ID
    inputBinding:
      prefix: --coords

  mirza:
    type: File
    doc: MIRZA output file
    inputBinding:
      prefix: --mirza

  contrafold:
    type: File
    doc: CONTRAfold output file
    inputBinding:
      prefix: --contrafold

  flanks:
    type: File
    doc: Flanks compositon output file
    inputBinding:
      prefix: --flanks

  distance:
    type: File
    doc: Distane to boundary output file
    inputBinding:
      prefix: --distance

  model_bls:
    type: File
    doc: Path to model with branch length score
    inputBinding:
      prefix: --model-bls

  model_nobls:
    type: File
    doc: Path to model without branch length score
    inputBinding:
      prefix: --model-nobls

  only_mirza:
    type:
    - "null"
    - type: enum
      symbols: ['yes', 'no']
    default: no
    doc: Calculate only MIRZA and DON'T calculate MIRZA BLS
    inputBinding:
      prefix: --only-mirza

  threshold:
    type: ["null", float]
    default: 0.12
    doc: Threshold for summing, defaults to 0.12
    inputBinding:
      prefix: --threshold

  split_by:
    type: ["null", string]
    default: NONE
    doc: If the header of fasta has multiple annotations eg. transcript_id|entrez_id|wikiname split it and take only one, defaults to NONE
    inputBinding:
      prefix: --split-by

  column:
    type: ["null", int]
    default: 0
    doc: 0 based column number to take after splitting, defaults to 0
    inputBinding:
      prefix: --column

  name:
    type: ["null", string]
    default: GeneID
    doc: Name of the id eg. gene, transcript etc , defaults to GeneID
    inputBinding:
      prefix: --name

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
      glob: $(inputs.coords.basename).score
