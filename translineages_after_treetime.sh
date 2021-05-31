#!/usr/bin/bash

tag="far20"
output_folder="../output/transmission_lineages/gennady/${tag}/treetime_noalign_median"
tree="${output_folder}/divergence_tree.nodates.newick"
#duplicates="../data/munched/gennady/merged_${tag}.duplicates" # merged_far20; by python merge_duplicates.py
nonrus_duplicates="../data/munched/gennady/definitely_all_nonrus_duplicates_of_nonrus_strains_${tag}" #  by python merge_nonrus_duplicates.py
meta="../data/raw/metadata_2021-03-12_09-12.tsv"
rpn_meta="../data/russian/meta_all.xlsx"

seqlength=29903

mkdir -p $output_folder

leaf_states_file=${output_folder}/leaf_states.csv
leaf_dates_file=${output_folder}/leaf_dates.csv
all_states_file=${output_folder}/Sankoff_states.csv
translin_file=${output_folder}/transmission_lineages.withduplicates.out
stripped_tree_file=${output_folder}/stripped_tree_with_labels_with_duplicates

# if [ -f "$tree" ]; then
    # echo "Using $tree that already exists"
# else 
    # echo "Trying to remove dates from nexus divergence tree to get $tree"
	# grep "Tree tree1" ${output_folder}/divergence_tree.nexus | tail -c +13 | head -c -2 >${output_folder}/divergence_tree.newick.temp
	# echo "(" >${output_folder}/divergence_tree.newick
	# cat ${output_folder}/divergence_tree.newick.temp >>${output_folder}/divergence_tree.newick
	# echo ");" >>${output_folder}/divergence_tree.newick
	# sed -r 's/\[&date=[0-9\.]+\]//g' ${output_folder}/divergence_tree.newick >$tree
# fi

# python meta_to_states.py --tree ${tree} --meta $meta --rpn_meta $rpn_meta --output $leaf_states_file # id	[2 if Russia 1 otherwise]	Country
# python meta_to_dates.py --tree ${tree} --meta $meta --rpn_meta $rpn_meta --output $leaf_dates_file

# echo "Estimating country states.."
# Rscript Sankoff_asr.R --tree ${tree} --states $leaf_states_file --output $all_states_file
# echo "Searching for russian transmission lineages.."
# python find_transmission_lineages.py --print_broader_subtree 5 --duplicates $nonrus_duplicates --tree ${all_states_file}.newick --states ${all_states_file}.probs  --countries $leaf_states_file --output $translin_file >${translin_file}.log
echo "Foreign strains timing.."
python -u timing.py --duplicates $nonrus_duplicates --tree ${all_states_file}.newick --dates $leaf_dates_file --countries $leaf_states_file  --entry_nodes_file ${translin_file}.entries --output $translin_file >$translin_file.date_test.log
# stats for dates inside transmission clusters
echo "Cluster timing.."
python cluster_timing.py --transmission_lineages_file $translin_file --dates $leaf_dates_file --output ${output_folder}/cluster_strains.dates_stats
python entry_to_date.py --tree ${all_states_file}.newick --tree_with_dates ${output_folder}/divergence_tree.nexus --output ${output_folder}/entry_to_date.csv