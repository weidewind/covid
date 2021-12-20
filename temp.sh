#!/usr/bin/bash
while getopts c: option
do
		case "${option}" in
			c) config_file=${OPTARG};;
		esac
done

source $config_file

set -e # exit if any command exits with non-zero status

masked_gisaid_msa=/export/home/popova/workspace/covid/data/raw/mmsa_2021-11-17/2021-11-17_masked.fa
variants_folder=/export/home/popova/workspace/covid/output/variant_analysis/${tag}
variants_file=${variants_folder}/taxon_mutations
source $conda_config


xvfb-run python -u strip_tree_on_closest.py --tree ${output_folder}/Sankoff_states.csv.newick --melted_foreign ${variants_file}_foreign \
--dates $leaf_dates_file --entry_nodes_file $variants_file --output ${variants_folder}/stripped_tree --scale 10


