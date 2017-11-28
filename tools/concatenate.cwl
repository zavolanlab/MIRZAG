#!/usr/bin/env cwl-runner

cwlVersion: v1.0

class: CommandLineTool
baseCommand: cat

inputs:

  files:
    type: File[]
    inputBinding:
      position: 1

  output_file:
    type: string
    default: concatenated_files

stdout: $(inputs.output_file)

outputs:
  output:
    type: stdout
