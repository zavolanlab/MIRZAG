#!/usr/bin/env cwl-runner

cwlVersion: v1.0

class: CommandLineTool
baseCommand: ['awk', '-f', 'script.awk']

requirements:
  InitialWorkDirRequirement:
    listing:
      - entryname: script.awk
        entry: |
          {print | "gzip > "$2".seedcount"}

inputs:
  file:
    type: File
    inputBinding:
      position: 1

outputs:
  output:
    type: File[]
    outputBinding:
      glob: "*.seedcount"
