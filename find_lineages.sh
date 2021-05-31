#!/usr/bin/bash

perl meta_to_states.pl --tree ../data/ft_SH.tree --meta ../data/metadata_2021-03-04_10-31.tsv --output ../data/ft_SH_states.csv
Rscript Sankoff_asr.R --tree ../data/ft_SH.tree --states ../data/ft_SH_states.csv --output ../data/ft_SH_Sankoff_output
python find_transmission_lineages.py --tree ../data/ft_SH_Sankoff_output.newick --states ../data/ft_SH_Sankoff_output.probs --output ../data/ft_SH_testpy >../data/ft_SH_pylog