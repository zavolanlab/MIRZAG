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

# TODO implement 'seed scan' (a parallel option for execution)

steps:

  prepare_mirnas:
    run: tools/prepare_mirnas.cwl
    in:
      input: input_mirna
    out:
      [ output ]

  generate_mrna_chunks:
    run: tools/generate_mrna_chunks.cwl
    in:
      input: input_mrna
    out:
      [ output ]

  mirna_expressions:
    run: tools/generate_mirna_expressions.cwl
    in:
      input: input_mirna
    out:
      [ output ]

  mirza:
    run: tools/mirza.cwl
    scatter: [ mrna, mirna ]
    scatterMethod: flat_crossproduct
    in:
      expressions: mirna_expressions/output
      mrna: generate_mrna_chunks/output
      mirna: prepare_mirnas/output
    out:
      [ output ]

  mirza_parser:
    run: tools/extract_data_from_mirza_output.cwl
    scatter: [ mirza ]
    in:
      mirza: mirza/output
      mrna: input_mrna
    out:
      [ output ]

  concatenate_mirza:
    run: tools/concatenate.cwl
    in:
      files: mirza_parser/output
      output_file:
        default: "mirza_scan_results"
    out:
      [ output ]

  filter_duplicates_from_scan:
    run: tools/filter_duplicates_from_scan.cwl
    in:
      coords: concatenate_mirza/output
    out:
      [ output ]

  split_mirza_results:
    run: tools/split_mirza_results.cwl
    in:
      file: filter_duplicates_from_scan/output
    out:
      [ output ]

  calculate_mirza:
    run: tools/calculate_MIRZA.cwl
    scatter: coords
    in:
      seq: input_mrna
      motifs: input_mirna
      coords: split_mirza_results/output
      tree: input_tree
      msa: input_multiple_alignments
    out:
      [ output ]

  calculate_contrafold:
    run: tools/calculate_contrafold.cwl
    scatter: coords
    in:
      seq: input_mrna
      coords: split_mirza_results/output
    out:
      [ output ]

  calculate_flanks:
    run: tools/calculate_flanks_composition.cwl
    scatter: coords
    in:
      seq: input_mrna
      coords: split_mirza_results/output
    out:
      [ output ]

  calculate_distance:
    run: tools/calculate_distance.cwl
    scatter: coords
    in:
      seq: input_mrna
      coords: split_mirza_results/output
    out:
      [ output ]

  merge_results_add_probability_and_calculate_per_gene_score:
    run: tools/merge_results_add_probability_and_calculate_per_gene_score.cwl
    scatter: [ coords, mirza, contrafold, flanks, distance ]
    scatterMethod: dotproduct
    in:
      coords: split_mirza_results/output
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
