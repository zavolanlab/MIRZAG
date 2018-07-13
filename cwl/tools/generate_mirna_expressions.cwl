#!/usr/bin/env cwl-runner

cwlVersion: v1.0

# TODO change Ruby to awk

class: CommandLineTool
baseCommand: ruby
hints:
  SoftwareRequirement:
    packages:
    - package: 'Ruby'
      version:
      - '2.3.3'
  DockerRequirement:
    dockerPull: ruby:2.4.2-slim

stdout: $(inputs.output_name)

arguments:
  - "-ne"
  - "puts \"#{$_.rstrip()[1..-1]}\t1\" if $_.start_with?(\">\")"

inputs:
  input:
    type: stdin
  output_name:
    type: string
    default: "expressions.tsv"

outputs:
  output:
    type: File
    outputBinding:
      glob: $(inputs.output_name)
