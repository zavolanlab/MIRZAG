#!/usr/bin/env cwl-runner

cwlVersion: v1.0

class: CommandLineTool
baseCommand: MIRZA

stdout: $(inputs.mirna.basename)_$(inputs.mrna.basename).mirza

arguments:
  - valueFrom: "50"
    position: 4
  - valueFrom: noupdate
    position: 5

inputs:

  expressions:
    type: File
    inputBinding:
      position: 1

  mrna:
    type: File
    inputBinding:
      position: 2

  mirna:
    type: File
    inputBinding:
      position: 3

outputs:

  output:
    type: File
    outputBinding:
      glob: $(inputs.mirna.basename)_$(inputs.mrna.basename).mirza
