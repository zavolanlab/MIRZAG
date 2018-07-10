#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow

requirements:
   - class: ScatterFeatureRequirement
   - class: SubworkflowFeatureRequirement

inputs:

  input_mirna: File
  input_mrna: File
  input_tree: File
  input_multiple_alignments: File
  input_model_with_bls: File
  input_model_without_bls: File
  output_file_name: string
  settings_split_by: string
  settings_index_after_split: int
  settings_mirza_threshold: float
  settings_contextLen_L: int
  settings_contextLen_U: int

outputs:

  mirzag_output:
    type: File
    outputSource: analyze/output

steps:

  prepare_mirnas:
    run: tools/prepare_mirnas.cwl
    in:
      input: input_mirna
    out:
      [ output ]

  split:
    run: tools/count_miRNA_seeds_and_filter_duplicates.cwl
    scatter: [ motifs ]
    in:
      motifs: prepare_mirnas/output
      seqs: input_mrna
      split_by: settings_split_by
      index_after_split: settings_index_after_split
    out:
      [ output ]

  analyze:
    run: tools/analysis.cwl
    in:
      input_mirna: input_mirna
      input_mrna: input_mrna
      input_coords:
        source: split/output
      input_tree: input_tree
      input_multiple_alignments: input_multiple_alignments
      input_model_with_bls: input_model_with_bls
      input_model_without_bls: input_model_without_bls
      output_file_name: output_file_name
      settings_split_by: settings_split_by
      settings_index_after_split: settings_index_after_split
      settings_mirza_threshold: settings_mirza_threshold
      settings_contextLen_L: settings_contextLen_L
      settings_contextLen_U: settings_contextLen_U
    out:
      [ output ]
