#!/usr/bin/bash

tree="../data/munched/gennady/for_iqtree.fasta.treefile"
duplicates="../data/munched/gennady/merged_far8.duplicates" # merged_far20; by python merge_duplicates.py
output_folder="../output/transmission_lineages/gennady/large"
meta="../data/raw/metadata_2021-03-12_09-12.tsv"
rpn_meta="../data/russian/meta_all.xlsx"



leaf_states_file=${output_folder}/leaf_states.csv
leaf_dates_file=${output_folder}/leaf_dates.csv
all_states_file=${output_folder}/Sankoff_states.csv
translin_file=${output_folder}/transmission_lineages.withduplicates.out
stripped_tree_file=${output_folder}/stripped_tree_with_labels_with_duplicates

python meta_to_states.py --tree ${tree}.rerooted --meta $meta --rpn_meta $rpn_meta --output $leaf_states_file 
python meta_to_dates.py --tree ${tree}.rerooted --meta $meta --rpn_meta $rpn_meta --output $leaf_dates_file
