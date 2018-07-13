#!/usr/bin/env cwl-runner

cwlVersion: v1.0

class: CommandLineTool
baseCommand: ['gzip', '-c']

inputs:
  file:
    type: File
    inputBinding:
      position: 1

stdout: $(inputs.file.basename).gz

outputs:
  output:
    type: File
    outputBinding:
      glob: $(inputs.file.basename).gz
