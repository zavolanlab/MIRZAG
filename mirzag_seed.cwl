#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow

requirements:
   - class: ScatterFeatureRequirement

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

outputs:

  mirzag_output:
    type: File
    outputSource: concatenate_results/output

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

  calculate_mirza:
    run: tools/calculate_MIRZA.cwl
    scatter: coords
    in:
      seq: input_mrna
      motifs: input_mirna
      coords: split/output
      tree: input_tree
      msa: input_multiple_alignments
    out:
      [ output ]

  calculate_contrafold:
    run: tools/calculate_contrafold.cwl
    scatter: coords
    in:
      seq: input_mrna
      coords: split/output
    out:
      [ output ]

  calculate_flanks:
    run: tools/calculate_flanks_composition.cwl
    scatter: coords
    in:
      seq: input_mrna
      coords: split/output
    out:
      [ output ]

  calculate_distance:
    run: tools/calculate_distance.cwl
    scatter: coords
    in:
      seq: input_mrna
      coords: split/output
    out:
      [ output ]

  merge_results_add_probability_and_calculate_per_gene_score:
    run: tools/merge_results_add_probability_and_calculate_per_gene_score.cwl
    scatter: [ coords, mirza, contrafold, flanks, distance ]
    scatterMethod: dotproduct
    in:
      coords: split/output
      mirza: calculate_mirza/output
      contrafold: calculate_contrafold/output
      flanks: calculate_flanks/output
      distance: calculate_distance/output
      model_bls: input_model_with_bls
      model_nobls: input_model_without_bls
      split_by: settings_split_by
      column: settings_index_after_split
    out:
      [ output ]

  concatenate_results:
    run: tools/concatenate.cwl
    in:
      files: merge_results_add_probability_and_calculate_per_gene_score/output
      output_file: output_file_name
    out:
      [ output ]
