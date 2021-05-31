#!/usr/bin/bash

tree="../data/munched/gennady/for_iqtree_far20.fasta.treefile"
meta="../data/raw/metadata_2021-03-12_09-12.tsv"
rpn_meta="../data/russian/meta_all.xlsx"
output_folder="../output/transmission_lineages/gennady"

mkdir -p $output_folder

leaf_states_file=${output_folder}/leaf_states.csv
all_states_file=${output_folder}/Sankoff_states.csv

python meta_to_states.py --tree $tree --meta $meta --rpn_meta $rpn_meta --output $leaf_states_file # id	[2 if Russia 1 otherwise]	Country