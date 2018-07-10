#!/usr/bin/env cwl-runner

cwlVersion: v1.0

class: CommandLineTool
baseCommand: rg_filter_duplicates_from_scan.py

doc: |

  Filter duplicates in coordinates by id, miRNA and sequence


inputs:

  coords:
    type: File
    doc: Coordinates file
    inputBinding:
      prefix: --coords

  split_by:
    type: string
    default: "NONE"
    doc: Split id by the string
    inputBinding:
      prefix: --split-by

  index_after_split:
    type: int
    default: 0
    doc: After split take this column as new id, 0 based
    inputBinding:
      prefix: --index-after-split

  verbose:
    type: ["null", boolean]
    default: False
    doc: Be loud!
    inputBinding:
      prefix: --verbose

stdout: $(inputs.coords.basename).filtered

outputs:
  output:
    type: File
    outputBinding:
      glob: $(inputs.coords.basename).filtered
