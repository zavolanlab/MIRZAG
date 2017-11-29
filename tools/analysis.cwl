#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow

requirements:
   - class: ScatterFeatureRequirement

inputs:

  input_mirna: File
  input_mrna: File
  input_coords: File[]
  input_tree: File
  input_multiple_alignments: File
  input_model_with_bls: File
  input_model_without_bls: File
  output_file_name: string
  settings_split_by: string
  settings_index_after_split: int

outputs:

  output:
    type: File
    outputSource: concatenate_results/output

steps:

  calculate_mirza:
    run: calculate_MIRZA.cwl
    scatter: coords
    in:
      seq: input_mrna
      motifs: input_mirna
      coords: input_coords
      tree: input_tree
      msa: input_multiple_alignments
    out:
      [ output ]

  calculate_contrafold:
    run: calculate_contrafold.cwl
    scatter: coords
    in:
      seq: input_mrna
      coords: input_coords
    out:
      [ output ]

  calculate_flanks:
    run: calculate_flanks_composition.cwl
    scatter: coords
    in:
      seq: input_mrna
      coords: input_coords
    out:
      [ output ]

  calculate_distance:
    run: calculate_distance.cwl
    scatter: coords
    in:
      seq: input_mrna
      coords: input_coords
    out:
      [ output ]

  merge_results_add_probability_and_calculate_per_gene_score:
    run: merge_results_add_probability_and_calculate_per_gene_score.cwl
    scatter: [ coords, mirza, contrafold, flanks, distance ]
    scatterMethod: dotproduct
    in:
      coords: input_coords
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
    run: concatenate.cwl
    in:
      files: merge_results_add_probability_and_calculate_per_gene_score/output
      output_file: output_file_name
    out:
      [ output ]
