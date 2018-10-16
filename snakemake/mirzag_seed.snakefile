rule count_miRNA_seeds_and_filter_duplicates:
    input:
        'results/split_mirnas/{part}.mirna.fa'
    output:
        'results/coords/{part}.seedcount'
    params:
        cluster_log = "results/cluster_log/count_miRNA_seeds_and_filter_duplicates_{part}.log"
    singularity:
        'docker://zavolab/mirzag-scripts:1'
    log:
        "results/local_log/count_miRNA_seeds_and_filter_duplicates_{part}.log"
    shell:
        "(rg_count_miRNA_seeds_and_filter_duplicates.py --motifs {input} \
        --seqs {config[input_mrna]} --split-by {config[settings_split_by]} \
        --index-after-split {config[settings_index_after_split]} \
        --output {output}) &> {log}"
