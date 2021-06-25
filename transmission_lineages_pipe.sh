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
translin_file=${output_folder}/transmission_lineages.withduplicates.out
stripped_tree_file=${output_folder}/stripped_tree_with_labels_with_duplicates
country_stat_file=${output_folder}/country_stat

mkdir -p $output_folder

# #python check_duplicates.py
# echo "Merging duplicates.."
# python merge_rus_duplicates.py --rtr_duplicates $rus_to_rus_duplicates --rtn_duplicates $rus_to_nonrus_duplicates --output $rus_duplicates
# python merge_nonrus_duplicates.py --nonrus_duplicates $nr_to_nr_duplicates --nonrus_clusters $nonrus_clusters --output $nonrus_duplicates
# python merge_nonrus_duplicates.py --nonrus_duplicates $nr_to_nr_duplicates --nonrus_clusters $nonrus_clusters --without_clusters --output $nonrus_duplicates_wo_clusters
# ##python check_merge.py # far20 hard-coded inside, change it to check other trees (and probably #todo make everything fall if it does not pass)

# echo "Rerooting and refluffing.."
# python meta_to_dates.py --tree ${tree} --output $leaf_dates_file
# python reroot.py --tree ${tree} --dates $leaf_dates_file --output ${tree}.rerooted
# python refluff.py --tree ${tree}.rerooted --duplicates $rus_duplicates --output ${tree}.rerooted.rus_refluffed

# python meta_to_states.py --tree ${tree}.rerooted.rus_refluffed --output $leaf_states_file # id	[2 if Russia 1 otherwise]	Country
# python meta_to_dates.py --tree ${tree}.rerooted.rus_refluffed --output $leaf_dates_file

# echo "State reconstruction.."
# Rscript Sankoff_asr.R --tree ${tree}.rerooted.rus_refluffed --states $leaf_states_file --output $all_states_file
# echo "Searching for transmission lineages.." 
#echo "python find_transmission_lineages.py --print_broader_subtree -1 --duplicates $nonrus_duplicates --tree ${all_states_file}.newick --states ${all_states_file}.probs  --countries $leaf_states_file --output $translin_file >${translin_file}.log"
# python find_transmission_lineages.py --print_broader_subtree 5 --duplicates $nonrus_duplicates --tree ${all_states_file}.newick --states ${all_states_file}.probs  --countries $leaf_states_file --output $translin_file >${translin_file}.log
#echo "Collecting country statistics.."
maxdist=0.00009
echo "python -u country_stat.py --duplicates $nonrus_duplicates --tree ${all_states_file}.newick --states ${all_states_file}.probs --countries $leaf_states_file --entry_nodes_file ${translin_file}.entries --max_distance ${maxdist} --output ${country_stat_file}.${maxdist} >${country_stat_file}.logger.${maxdist}"
#python -u country_stat.py --duplicates $nonrus_duplicates --tree ${all_states_file}.newick --states ${all_states_file}.probs --countries $leaf_states_file --entry_nodes_file ${translin_file}.entries --max_distance ${maxdist} --output ${country_stat_file}.${maxdist} >${country_stat_file}.logger.${maxdist}

#echo "Stripping tree.." 
#echo "xvfb-run python strip_tree_with_maxdepth.py --duplicates $nonrus_duplicates --tree ${all_states_file}.newick --states ${all_states_file}.probs --countries $leaf_states_file --entry_nodes_file ${translin_file}.entries --countries $leaf_states_file --max_distance 0.0001 --output $stripped_tree_file >${stripped_tree_file}.logger"
#xvfb-run python strip_tree_with_maxdepth.py --duplicates $nonrus_duplicates --tree ${all_states_file}.newick --states ${all_states_file}.probs --countries $leaf_states_file --entry_nodes_file ${translin_file}.entries --countries $leaf_states_file --max_distance 0.0001 --output $stripped_tree_file >${stripped_tree_file}.logger


# ## Internal nodes dating
# # trying to get dates of foreign and russian strains, which are closest to entry node (entry node dating for the poor)
#echo "Foreign strains timing.."
#echo "python -u timing.py --duplicates $nonrus_duplicates_wo_clusters --tree ${all_states_file}.newick --dates $leaf_dates_file --countries $leaf_states_file  --entry_nodes_file ${translin_file}.entries --output $translin_file >$translin_file.date_test.log"
#python -u timing.py --duplicates $nonrus_duplicates_wo_clusters --tree ${all_states_file}.newick --dates $leaf_dates_file --countries $leaf_states_file  --entry_nodes_file ${translin_file}.entries --output $translin_file >$translin_file.date_test.log
# # stats for dates inside transmission clusters
 #echo "Cluster timing.."
 #echo " python cluster_timing.py --transmission_lineages_file $translin_file --dates $leaf_dates_file --output ${output_folder}/cluster_strains.dates_stats"
 #python cluster_timing.py --transmission_lineages_file $translin_file --dates $leaf_dates_file --output ${output_folder}/cluster_strains.dates_stats


# # Treetime needs exactly one date for a leaf. Thus, for strains with duplicates, print the earliest/the latest/median defined date
# echo "Reprinting dates from meta to leaf.dates files.."
# python meta_to_dates_with_duplicates.py --tree ${all_states_file}.newick --duplicates $nonrus_duplicates --type min --output ${output_folder}/leaf_dates.csv.cleaned
# python meta_to_dates_with_duplicates.py --tree ${all_states_file}.newick --duplicates $nonrus_duplicates --type median --output ${output_folder}/leaf_dates.csv.cleaned
# python meta_to_dates_with_duplicates.py --tree ${all_states_file}.newick --duplicates $nonrus_duplicates --type all --output ${output_folder}/leaf_dates.csv.cleaned

# python meta_to_dates_with_duplicates.py --tree ${all_states_file}.newick  --duplicates $nonrus_duplicates_wo_clusters --type min --output ${output_folder}/leaf_dates.csv.cleaned.wo_clusters
# python meta_to_dates_with_duplicates.py --tree ${all_states_file}.newick  --duplicates $nonrus_duplicates_wo_clusters --type median --output ${output_folder}/leaf_dates.csv.cleaned.wo_clusters
# python meta_to_dates_with_duplicates.py --tree ${all_states_file}.newick  --duplicates $nonrus_duplicates_wo_clusters --type all --output ${output_folder}/leaf_dates.csv.cleaned.wo_clusters

# could not make treetime work under screen, hence nohup (does not work from script):
#echo "nohup unbuffer treetime --keep-root --sequence-length ${seqlength} --tree ${all_states_file}.newick --dates ${output_folder}/leaf_dates.csv.cleaned.wo_clusters.min  --outdir ${output_folder}/treetime_noalign_min >${output_folder}/treetime_noalign_min/treetime.logs &"
##nohup unbuffer treetime --keep-root --sequence-length $seqlength --tree ${all_states_file}.newick --dates ${output_folder}/leaf_dates.csv.cleaned.min  --outdir ${output_folder}/treetime_noalign_min >${output_folder}/treetime_noalign_min/treetime.logs &
##nohup unbuffer treetime --keep-root --sequence-length $seqlength --tree ${all_states_file}.newick --dates ${output_folder}/leaf_dates.csv.cleaned.min  --outdir ${output_folder}/treetime_noalign_median >${output_folder}/treetime_noalign_median/treetime.logs &
