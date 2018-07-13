# FIXME Snakemake has problems with double wild cards in dynamic(), hence
# this rule is only parsing the file, it doesn't split it into multiple outputs
rule generate_mrna_chunks:
    input:
        config['input_mrna']
    output:
        'results/part_1.mrna.fa'
    singularity:
        'docker:/zavolab/mirzag-scripts:1'
    shell:
        'rg_generate_utr_chunks.py --input {input} '
        '--output-dir results --part-size 90000000000'

rule mirna_expressions:
    input:
        config['input_mirna']
    output:
        'results/mirna_expressions.tsv'
    run:
        with open(input[0]) as fin, open(output[0], 'w') as fout:
            for line in fin:
                if line.startswith('>'):
                    fout.write('{}\t1\n'.format(line.rstrip()[1:]))

rule mirza:
    input:
        expressions = 'results/mirna_expressions.tsv',
        mrna = 'results/part_1.mrna.fa',
        mirna = 'results/split_mirnas/{part}.mirna.fa'
    output:
        'results/mirza_scan/{part}.mirza'
    singularity:
        'docker://zavolab/mirza:1'
    shell:
        'MIRZA {input} 50 noupdate > {output}'

rule mirza_parser:
    input:
        mirza = 'results/mirza_scan/{part}.mirza',
        seqs = config['input_mrna']
    output:
        'results/mirza_scan_parsed/{part}.tsv',
    singularity:
        'docker://zavolab/mirzag-scripts:1'
    shell:
        'rg_extract_data_from_mirza_output.py --input {input.mirza} '
        '--seqs {input.seqs} --threshold 50 --context 50 '
        '--output {output}'

rule concatenate_mirza:
    input:
        dynamic('results/mirza_scan_parsed/{part}.tsv')
    output:
        'results/mirza_scan_results.tsv'
    shell:
        'cat {input} > {output}'

rule filter_duplicates_from_scan:
    input:
        'results/mirza_scan_results.tsv'
    output:
        'results/mirza_scan_results.filtered.tsv'
    singularity:
        'docker://zavolab/mirzag-scripts:1'
    shell:
        'rg_filter_duplicates_from_scan.py --coords {input} '
        '--split-by NONE --index-after-split 0 > {output}'

rule split_mirza_results:
    input:
        'results/mirza_scan_results.filtered.tsv'
    output:
        dynamic('results/coords/{part}.seedcount')
    shell:
        "awk '{{print | \"gzip > results/coords/\"$2\".seedcount\"}}' {input}"
