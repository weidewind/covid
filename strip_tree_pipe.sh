#!/usr/bin/bash

tag="may31"
tree="../data/munched/gennady/${tag}/21_05_2021_for_iqtree_3_-6.fasta.treefile"
rus_to_rus_duplicates="../data/munched/gennady/${tag}/21_05_2021_rus_to_rus_dups.tsv" 
rus_to_nonrus_duplicates="../data/munched/gennady/${tag}/21_05_2021_rus_to_nonrus_dups.tsv"
rus_duplicates="../data/munched/gennady/${tag}/21_05_2021_all_rus_duplicates" 
nr_to_nr_duplicates="../data/munched/gennady/${tag}/21_05_2021_nonrus_to_nonrus_dups.tsv" #  by python merge_nonrus_duplicates.py
nonrus_clusters="../data/munched/gennady/${tag}/21_05_2021_nonrus_clusters_3_-6.list"
nonrus_duplicates="../data/munched/gennady/${tag}/21_05_2021_all_nonrus_duplicates"
# nonrus_duplicates_wo_clusters is used in timing.py
nonrus_duplicates_wo_clusters="../data/munched/gennady/${tag}/21_05_2021_nonrus_duplicates_wo_clusters"
output_folder="../output/transmission_lineages/gennady/${tag}"
seqlength=29903

leaf_states_file=${output_folder}/leaf_states.csv
leaf_dates_file=${output_folder}/leaf_dates.csv
all_states_file=${output_folder}/Sankoff_states.csv
translin_file=${output_folder}/transmission_lineages.withduplicates.all.out
nonrus_duplicates_after_defluffing="../data/munched/gennady/${tag}/21_05_2021_all_nonrus_duplicates_after_defluffing"

## add strains that were removed by defluffing (one foreign duplicate of russian strain is preserved in the tree, and all others are now added to duplicates file as duplicates of the preserved strain )

#python defluff.py --tree ${all_states_file}.newick --duplicates $rus_to_nonrus_duplicates --output ${output_folder}/defluffed
#cp $nonrus_duplicates $nonrus_duplicates_after_defluffing
#cat "${output_folder}/defluffed.duplicates" >>$nonrus_duplicates_after_defluffing
xvfb-run python strip_tree_on_closest.py --tree ${output_folder}/defluffed.newick --dates $leaf_dates_file --countries $leaf_states_file --duplicates $nonrus_duplicates_after_defluffing --entry_nodes_file ${translin_file}.entries --output ${output_folder}/all.defluffed_and_stripped
