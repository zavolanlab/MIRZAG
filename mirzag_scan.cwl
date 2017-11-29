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

  split:
    run: tools/split_mirza_results.cwl
    in:
      file: filter_duplicates_from_scan/output
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
    out:
      [ output ]
