Pipeline flow
*************

1. Prepare miRNAs for pipeline
==============================
.. argparse::
    :ref: scripts.rg_prepare_mirnas_for_mirza_and_split.parser
    :prog: rg_prepare_mirnas_for_mirza_and_split

2. Chunk UTRs
=============
.. argparse::
    :ref: scripts.rg_generate_utr_chunks.parser
    :prog: rg_generate_utr_chunks

3. Generate miRNA expressions
=============================
Generate expressions of miRNAs. This file is required by MIRZA and is composed just from the
miRNA ID and its expression. Here we set up expression for each miRNA to 1.

There is no special script dedicated to this function. It is just the bash command:

.. code-block:: bash

    cat input | ruby -ne 'puts "#{$_.rstrip()[1..-1]}\t1" if $_.start_with?(">")' > output


4a. Run MIRZA analysis
=======================
If the option "scan" as a miRNA target search is chosen putative targets are selected based on MIRZA interaction
energy. Thus, this is the  part of the pipeline scans the 3'UTRs with MIRZA to find it.

.. argparse::
    :ref: run_mirza_scan.parser
    :prog: run_mirza_scan

I. Calculate coordinates with MIRZA
-----------------------------------

Here the MIRZA algorithm is used to calculate the energy between miRNA and the 3'UTR fragments generated
before. MIRZA is launched from the command:

.. code-block:: bash

    MIRZA expressions.tab mrnas.fa mirnas.fa 50 noupdate

And the result is piped to the analysis script:

.. argparse::
    :ref: scripts.rg_extract_data_from_mirza_output.parser
    :prog: rg_extract_data_from_mirza_output

II. Merge and filter results
----------------------------

The results are merged with bash command and duplicate results resulting from eg.
overlapping 3'UTR fragments or transcripts from the same gene.

.. argparse::
    :ref: scripts.rg_filter_duplicates_from_scan.parser
    :prog: rg_filter_duplicates_from_scan

4b. Run seed scan
=================

This part is launched if the "seed" options is chosen when starting the pipeline. The putative
miRNA targets are found by simple seed match.

.. argparse::
    :ref: scripts.rg_count_miRNA_seeds_and_filter_duplicates.parser
    :prog: rg_count_miRNA_seeds_and_filter_duplicates


5. Run features analysis
========================

This is the main part of the pipeline where all the features are calculated.

I. Calculate MIRZA
------------------

.. argparse::
    :ref: scripts.rg_calculate_MIRZA.parser
    :prog: rg_calculate_MIRZA

II. Calculate accessibility with CONTRAfold
-------------------------------------------

.. argparse::
    :ref: scripts.rg_calculate_contrafold.parser
    :prog: rg_calculate_contrafold


III. Calculate flanks composition
---------------------------------

.. argparse::
    :ref: scripts.rg_calculate_flanks_composition.parser
    :prog: rg_calculate_flanks_composition

IV. Calculate distance to the boundary
--------------------------------------

.. argparse::
    :ref: scripts.rg_calculate_distance.parser
    :prog: rg_calculate_distance


6. Merge and add probabilities
==============================

.. argparse::
    :ref: scripts.rg_merge_results_add_probability_and_calculate_per_gene_score.parser
    :prog: rg_merge_results_add_probability_and_calculate_per_gene_score


7. Collect results
==================

In the end all the results are collected with bash command:

.. code-block:: bash

    zcat output_dir/\*.score > cwd/mirza_g_results_protocol.tab

