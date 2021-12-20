#!/usr/bin/bash
while getopts c: option
do
		case "${option}" in
			c) config_file=${OPTARG};;
		esac
done


source $config_file

window=14
melted_entries=${output_folder}/melted_entries.withdates.csv

mkdir -p $output_folder

set -e # exit if any command exits with non-zero status

# Rscript melt_entries.R  --pangolined_file $rus_withdups_pangolined  \
  # --entry_strains_file ${translin_file} --leaf_dates_file $leaf_dates_file \
  # --output $melted_entries --dates_stats_file ${translin_file}.dates_stats --entry_dates_file ${translin_file}.dates_stats.corrected.csv

Rscript proliferation_multilin.R  --window $window --melted_entries_file $melted_entries \
--important_lineages $important_lineages --output ${output_folder}/strains_to_entry