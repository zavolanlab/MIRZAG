Usage
*****

Basic usage
===========

Command to launch the pipeline is as follows:

.. code-block:: bash

    python MIRZA_G_pipeline.py run --config congig.ini --name-suffix name_of_the_run

All parameters for the script:

.. argparse::
    :ref: MIRZA_G_pipeline.parser
    :prog: MIRZA-G


Preparing config file
=====================

Copy config_example.ini from MIRZA-G directory to your working directory (directory
where you want to perform calculation, WD):

.. code-block:: bash

    cd Your/Working/Direcory
    cp Path/To/MIRZA-G/config_example.ini config.ini

Set all the necessary paths in your config.ini file as indicated in the comments inside the file. The most importand are:
 * **motifs**: "Path/To/miRNAs.fa" - abs path to an input fasta file with mi/siRNA sequences of length 21 or more
 * **seqs**: "Path/To/MIRZA-G/data/UTR_Sequences.fa" - abs path to a fasta file with the UTR sequences from which the coordinate file will be generated (you can use 3'UTR sequences in the
   MIRZA-G/data directory, for this file there are also alignments for conservation precalculated)
 * **mirza_binary**: "MIRZA" - path to MIRZA binary (or how you invoke it in the bash)
 * **contrafold_binary**: "contrafold" - path to CONTRAfold binary (or how you invoke in the bash)

Models paths:
 * **model_with_bls**: "Path/To/MIRZA-G/data/glm-with-bls.bin" - abs path to the model with BLS (you can find it in the pipeline/data directory)
 * **model_without_bls**: "Path/To/MIRZA-G/data/glm-without-bls.bin" - same as before

Additionally when you would like to calculate with evolutionary conservation you have to make sure that the variable run_only_MIRZA in CalculateMIRZA task is set to “no” instead of “yes” and that you provide proper paths with aligned UTRs and evolutionary tree:
 * **phylogenetic_tree**: "Path/To/MIRZA-G/data/human_tree.nh" - abspath to provided phylogenetic tree
 * **alignment_directory**: "Path/To/MIRZA-G/data/HumanAlignments/" - abspath to provided human alignments directory. If you downloaded package from CLIPz website
      this directory is already in the MIRZA-G directory. If you downloaded from GitHub you have to download it additionally.

If you would like to run it on cluster follow instructions in the configuration file and ask your admin what parameters you need to set
up before (like DRMAA path, modules necessary, queues names etc.). All these parameters can be set up in config.ini.

To run it locally it takes ~70 to 90 seconds for one miRNA without conservation calculation and ~170 seconds with calculation (This
might be substantial amount of time (up to half an hour per miRNA) for worse processors).


Example
=======

To test the pipeline go to the tests directory and run:

.. code-block:: bash

    cd Path/To/MIRZA-G/tests
    bash rg_run_test.sh help

.. note::

    Usage: rg_run_test.sh clean/run [MIRZA/binary/path] ['CONTRAfold/binary/path']

And if you have installed MIRZA and CONTRAfold to default locations (MIRZA and contrafold) run:

.. code-block:: bash

    bash rg_run_test.sh run

Otherwise provide paths to **BOTH** of them:

.. code-block:: bash

    bash rg_run_test.sh run Path/To/MIRZA/binary Path/To/CONTRAfold/binary

