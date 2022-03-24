#!/usr/bin/bash
while getopts c: option
do
		case "${option}" in
			c) config_file=${OPTARG};;
		esac
done

source $config_file
source timing.sh
source closest.sh

mkdir -p $output_folder
cp $config_file ${output_folder}/$config_file

set -e # exit if any command exits with non-zero status



## $1 is imports_dates_stats
## $2 is exports_dates_stats
## $3 output file name
## $4 reverse: 1 if for exports file we need to change rus and foreign, 0 otherwise
function compare {
	echo ${output_folder}/$3
	Rscript compare_imports_to_exports.R --imports_dates_stats $1 \
--exports_dates_stats $2 --reverse $4 \
--entries_file ${translin_file} --output_file ${output_folder}/$3 \
--imports_data ${output_folder}/${tag}_data.csv \
--lineages B.1.1.523,B.1.1.317,AY.122,AT.1,B.1.1.397,B.1.1.7,BA

	Rscript compare_imports_to_exports.R --imports_dates_stats $1 \
--exports_dates_stats $2 --reverse $4 \
--entries_file ${translin_file} --output_file ${output_folder}/$3
}

translin_file=${translin_file}.broken

# entries_ds=$(closest ${translin_file}.entries)
# exports_ds=$(closest ${translin_file}.exports)
# compare $entries_ds $exports_ds imports_to_exports.inout 1

random_entries_ds=$(closest_for_random ${output_folder}/withleaves_middle_entries_ 3000)
random_exports_ds=$(closest_for_random ${output_folder}/withleaves_middle_exports_ 3000)
compare $random_entries_ds $random_exports_ds withleaves_middle_random_imports_to_exports.inout 1



# echo "compare_imports_to_exports.."
# #compare tree_distance_to_closest_foreign_strain for exports and primary and secondary imports 
# Rscript compare_imports_to_exports.R --imports_dates_stats ${translin_file}.entries.dates_stats \
# --exports_dates_stats ${translin_file}.exports.dates_stats \
# --entries_file ${translin_file} --output_file ${output_folder}/imports_to_exports \
# --imports_data ${output_folder}/${tag}_data.csv \
# --lineages B.1.1.523,B.1.1.317,AY.122,AT.1,B.1.1.397,B.1.1.7,BA

# Rscript compare_imports_to_exports.R --imports_dates_stats ${translin_file}.dates_stats \
# --exports_dates_stats ${translin_file}.exports.dates_stats \
# --entries_file ${translin_file} --output_file ${output_folder}/imports_to_exports

# echo "country_stat.py.."
# parallel python country_stat.py --tree ${all_states_file}.newick.broken --countries $leaf_states_file --max_distance {} \
# --entry_nodes_file ${translin_file}.entries --include_inner --output ${output_folder}/country_stat.include_inner.{} ::: 0 1 2 3

# parallel python country_stat.py --tree ${all_states_file}.newick.broken --countries $leaf_states_file --max_distance {} \
# --entry_nodes_file ${translin_file}.exports --include_inner --output ${output_folder}/country_stat.include_inner.exports.{} ::: 0 1 2 3

# echo "cluster_nodes_withdist.py.."
# python cluster_nodes_withdist.py --tree ${all_states_file}.newick.broken --entries_file ${translin_file}.entries \
# --exports_file ${translin_file}.entries.exports --output ${output_folder}/translin.inner.withdist
# grep "internal" ${output_folder}/translin.inner.withdist | cut -d, -f2 >${output_folder}/translin.inner
# parallel python country_stat.py --tree ${all_states_file}.newick.broken --countries $leaf_states_file --max_distance {} \
# --entry_nodes_file ${output_folder}/translin.inner --include_inner --output ${output_folder}/country_stat.include_inner.inner.{} ::: 0 1 2 3

# echo "entry_parents_withdist.py.."
# python entry_parents_withdist.py --tree ${all_states_file}.newick.broken --steps 1 \
# --entries_file ${translin_file}.entries --output ${output_folder}/translin.parents.1.withdist
# tail -n +2 ${output_folder}/translin.parents.1.withdist | cut -d, -f2  >${output_folder}/translin.parents.1
# parallel python country_stat.py --tree ${all_states_file}.newick.broken --countries $leaf_states_file --max_distance {} \
# --entry_nodes_file ${output_folder}/translin.parents.1 --include_inner --output ${output_folder}/country_stat.include_inner.parents.1.{} ::: 0 1 2 3

# echo "rus_to_foreign_stats.R.."
# Rscript rus_to_foreign_stats.R --folder ${output_folder}/ --output ${output_folder}/compare



## compare distances from entry to closest ingroup and closest outgroup 




