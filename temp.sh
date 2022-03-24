#!/usr/bin/bash
while getopts c: option
do
		case "${option}" in
			c) config_file=${OPTARG};;
		esac
done

source $config_file
source timing.sh

mkdir -p $output_folder
cp $config_file ${output_folder}/$config_file

set -e # exit if any command exits with non-zero status

translin_file=${translin_file}.broken

echo "compare_imports_to_exports.."
# compare tree_distance_to_closest_foreign_strain for exports and primary and secondary imports 
Rscript compare_imports_to_exports.R --imports_dates_stats ${translin_file}.dates_stats \
--exports_dates_stats ${translin_file}.exports.dates_stats \
--entries_file ${translin_file} --output_file ${output_folder}/imports_to_exports \
--imports_data ${output_folder}/${tag}_data.csv --lineages B.1.1.523,B.1.1.317,AY.122,B.1.1.397,B.1.1.7