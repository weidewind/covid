#!/usr/bin/bash

tree="../data/munched/gennady/for_iqtree.fasta.treefile"
duplicates="../data/munched/gennady/merged_far8.duplicates" # merged_far20; by python merge_duplicates.py
meta="../data/raw/metadata_2021-03-12_09-12.tsv"
rpn_meta="../data/russian/meta_all.xlsx"
output_folder="../output/transmission_lineages/gennady/large"
seqlength=29903

mkdir -p $output_folder

leaf_states_file=${output_folder}/leaf_states.csv
leaf_dates_file=${output_folder}/leaf_dates.csv
all_states_file=${output_folder}/Sankoff_states.csv
translin_file=${output_folder}/transmission_lineages.withduplicates.out
stripped_tree_file=${output_folder}/stripped_tree_with_labels_with_duplicates


#python meta_to_states.py --tree ../data/munched/gennady/for_iqtree.fasta.treefile.rerooted --meta ../data/raw/metadata_2021-03-12_09-12.tsv --rpn_meta ../data/russian/meta_all.xlsx --output ../output/transmission_lineages/gennady/large/leaf_states.csv
python meta_to_states.py --tree ${tree} --meta $meta --rpn_meta $rpn_meta --output $leaf_states_file # id	[2 if Russia 1 otherwise]	Country
python meta_to_dates.py --tree ${tree} --meta $meta --rpn_meta $rpn_meta --output $leaf_dates_file

python reroot.py --tree $tree --output ${tree}.rerooted

Rscript Sankoff_asr.R --tree ${tree}.rerooted --states $leaf_states_file --output $all_states_file
#python find_transmission_lineages.py --tree ${all_states_file}.newick --states ${all_states_file}.probs  --output ${output_folder}/transmission_lineages.out >${output_folder}/transmission_lineages.log
#python find_transmission_lineages.py --duplicates /export/home/popova/workspace/covid/data/munched/gennady/merged_far8.duplicates --tree ../output/transmission_lineages/gennady/large/Sankoff_states.csv.newick --states ../output/transmission_lineages/gennady/large/Sankoff_states.csv.probs  --output ../output/transmission_lineages/gennady/large/transmission_lineages.withduplicates.out >../output/transmission_lineages/gennady/large/transmission_lineages.withduplicates.log
python find_transmission_lineages.py --duplicates $duplicates --tree ${all_states_file}.newick --states ${all_states_file}.probs  --countries $leaf_states_file --output $translin_file >${translin_file}.log
xvfb-run python strip_tree_with_maxdepth.py --duplicates $duplicates --tree ${all_states_file}.newick --states ${all_states_file}.probs --countries $leaf_states_file --entry_nodes_file ${translin_file}.entries --countries $leaf_states_file --max_distance 0.0001 --output $stripped_tree_file >${stripped_tree_file}.logger

## Internal nodes dating
# trying to get dates of foreign and russian strains, which are closest to entry node (entry node dating for the poor)
python -u timing.py --duplicates $duplicates --tree ${all_states_file}.newick --dates $leaf_dates_file --countries $leaf_states_file  --entry_nodes_file ${translin_file}.entries --output $translin_file >$translin_file.date_test.log
# stats for dates inside transmission clusters
python cluster_timing.py --transmission_lineages_file ../output/transmission_lineages/gennady/far20/transmission_lineages.withduplicates.out --dates ../output/transmission_lineages/gennady/far20/leaf_dates.csv --output ../output/transmission_lineages/gennady/far20/cluster_strains.dates_stats

# for strains with duplicates, print the earliest defined date
python meta_to_dates_with_duplicates.py --tree ../data/munched/gennady/for_iqtree_far20.fasta.treefile.rerooted  --meta ../data/raw/metadata_2021-03-12_09-12.tsv --rpn_meta ../data/russian/meta_all.xlsx --duplicates ../data/munched/gennady/merged_far20.duplicates --type min --output ../output/transmission_lineages/gennady/far20/leaf_dates.csv.cleaned
# could not make treetime work under screen, hence nohup:
nohup unbuffer treetime --tree /export/home/popova/workspace/covid/data/munched/gennady/for_iqtree_far20.fasta.treefile --dates /export/home/popova/workspace/covid/output/transmission_lineages/gennady/far20/leaf_dates.csv.cleaned.min --aln  /export/home/popova/workspace/covid/data/munched/gennady/for_iqtree_far20.fasta --outdir /export/home/popova/workspace/covid/output/transmission_lineages/gennady/far20/treetime_min >/export/home/popova/workspace/covid/output/transmission_lineages/gennady/far20/treetime_min/treetime.logs &
# for treetime_all (consider all dates, not only mean or min) we need to add all duplicates into the tree:
python refluff.py --tree ../data/munched/gennady/for_iqtree_far20.fasta.treefile.rerooted --duplicates ../data/munched/gennady/merged_far20.duplicates --output ../data/munched/gennady/for_iqtree_far20.fasta.treefile.rerooted.refluffed

nohup unbuffer treetime --keep-root --sequence-length $seqlength --tree /export/home/popova/workspace/covid/data/munched/gennady/for_iqtree_far20.fasta.treefile.rerooted.refluffed --dates /export/home/popova/workspace/covid/output/transmission_lineages/gennady/far20/leaf_dates.csv.cleaned.all  --outdir /export/home/popova/workspace/covid/output/transmission_lineages/gennady/far20/treetime_noalign_all >/export/home/popova/workspace/covid/output/transmission_lineages/gennady/far20/treetime_noalign_all/treetime.logs &

screen -U -S treedater_far20 bash -c 'Rscript run_treedater.R --tree ../data/munched/gennady/for_iqtree_far20.fasta.treefile.rerooted --dates /export/home/popova/workspace/covid/output/transmission_lineages/gennady/far20/leaf_dates.csv.mean --upper /export/home/popova/workspace/covid/output/transmission_lineages/gennady/far20/leaf_dates.csv.max --lower /export/home/popova/workspace/covid/output/transmission_lineages/gennady/far20/leaf_dates.csv.min --output /export/home/popova/workspace/covid/output/transmission_lineages/gennady/far20/treedater --seqlength $seqlength >/export/home/popova/workspace/covid/output/transmission_lineages/gennady/far20/treedater/logs_and_summary'