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
deprecated_masked_gisaid_msa=/export/home/popova/workspace/covid/data/raw/mmsa_2021-09-23/mmsa_2021-09-23/2021-09-23_masked.fa
variants_folder=/export/home/popova/workspace/covid/output/variant_analysis/${tag}
variants_file=${variants_folder}/taxon_mutations
entry_mut_stats=${variants_folder}/entry_mut_stats
source $conda_config

mkdir -p $variants_folder

## $rus_withdupls_pangolined are created by rus_withdupls_pangolined.sh

# conda activate sarscov2phylo
# tail -n +2 $rus_withdupls_pangolined | cut -d, -f1 >${output_folder}/rus_ids.list
# faFilter -namePatList=${output_folder}/rus_ids.list $masked_gisaid_msa ${output_folder}/rus.gisaid.mmsa
# conda deactivate

# python get_strains_with_mutations.py --fasta ${output_folder}/rus.gisaid.mmsa,/export/home/popova/workspace/covid/data/raw/raw_rus_fasta_30sep/2021-09-30_all_rus_for_usher_nextstrain.fasta \
# --pangolined $rus_withdupls_pangolined \
# --melted_entries ${output_folder}/melted_entries.withdates.csv \
# --output $variants_file \
# --tree $tree \
# --dates ${output_folder}/leaf_dates.csv

# Rscript get_melted_foreign_strains_for_entries.R --dates_stats ${output_folder}/transmission_lineages.withduplicates.out.dates_stats \
# --variants_file $variants_file --output ${variants_folder}/melted_foreign

## find ids missing in the previous version of ${output_folder}/close_foreign.pangolined
# tail -n +2 ${variants_folder}/melted_foreign | cut -d';' -f3 >${output_folder}/all_close_foreign_ids.list
# sort ${output_folder}/all_close_foreign_ids.list | uniq >${output_folder}/uniq_close_foreign_ids.list
# >${output_folder}/close_foreign_ids.list
# while IFS="" read -r id || [ -n "$id" ] # otherwise it will skip the trailing line
# do
  # grep "${id}," ${output_folder}/close_foreign.pangolined || echo $id >>${output_folder}/close_foreign_ids.list
# done <${output_folder}/uniq_close_foreign_ids.list



# conda activate sarscov2phylo
# faFilter -namePatList=${output_folder}/close_foreign_ids.list $masked_gisaid_msa ${output_folder}/close_foreign_additional.gisaid.mmsa
# ##find ids missing in the new gisaid mmsa, look for them in the previous version of mmsa
# >${output_folder}/deprecated_close_foreign_ids.list
# while IFS="" read -r id || [ -n "$id" ] # otherwise it will skip the trailing line
# do
  # grep "${id}" ${output_folder}/close_foreign_additional.gisaid.mmsa || echo $id >>${output_folder}/deprecated_close_foreign_ids.list
# done <${output_folder}/close_foreign_ids.list
# faFilter -namePatList=${output_folder}/deprecated_close_foreign_ids.list $deprecated_masked_gisaid_msa ${output_folder}/close_foreign_additional_deprecated.gisaid.mmsa
# cat ${output_folder}/close_foreign_additional_deprecated.gisaid.mmsa >>${output_folder}/close_foreign_additional.gisaid.mmsa
# conda deactivate
# conda activate pangolin
# pangolin ${output_folder}/close_foreign_additional.gisaid.mmsa --outfile ${output_folder}/close_foreign_additional.pangolined
# conda deactivate
# tail -n +2 ${output_folder}/close_foreign_additional.pangolined >>${output_folder}/close_foreign.pangolined

# cat ${output_folder}/close_foreign_additional.gisaid.mmsa >>${output_folder}/close_foreign.gisaid.mmsa

# python get_foreign_strains_with_mutations.py --fasta ${output_folder}/close_foreign.gisaid.mmsa \
# --melted_entries ${variants_folder}/melted_foreign \
# --pangolined ${output_folder}/close_foreign.pangolined \
# --tree $tree \
# --output ${variants_file}_foreign 

# python entry_mut_stats.py --variants_file $variants_file --output ${variants_folder}/entry_mut_stats

# thr=0.33

# Rscript merge_data_for_plots.R --entry_weight_file ${output_folder}/entry_freq.coverage_${thr}.csv \
# --corrected_dates_file ${translin_file}.dates_stats.corrected.csv \
# --entry_dates_file ${translin_file}.dates_stats --pangolined_file $rus_withdupls_pangolined \
# --cluster_dates_file ${output_folder}/cluster_strains.dates_stats --entry_strains_file ${translin_file} \
# --important_lineages $important_lineages --output ${output_folder}/${tag}_data.csv --coverage $thr

# Rscript rug_plots.R --data_file ${output_folder}/${tag}_data.csv --mutations_file $variants_file --foreign_mutations_file ${variants_file}_foreign \
# --melted_entries_file ${output_folder}/melted_entries.withdates.csv --districts ${output_folder}/leaf_districts.csv \
# --output_folder $output_folder --tag $tag --coverage $thr

# parallel python country_stat.py --tree ${output_folder}/Sankoff_states.csv.newick --countries $leaf_states_file --max_distance {} \
# --entry_nodes_file ${translin_file}.entries --include_inner --output ${output_folder}/country_stat.include_inner.{} ::: 0 1 2 3

# parallel python country_stat.py --tree ${output_folder}/Sankoff_states.csv.newick --countries $leaf_states_file --max_distance {} \
# --entry_nodes_file ${output_folder}/translin.exports --include_inner --output ${output_folder}/country_stat.include_inner.exports.{} ::: 0 1 2 3


# python cluster_nodes_withdist.py --tree ${output_folder}/Sankoff_states.csv.newick --entries_file ${translin_file}.entries \
# --exports_file ${output_folder}/translin.entries.exports --output ${output_folder}/translin.inner.withdist

# grep "internal" ${output_folder}/translin.inner.withdist | cut -d, -f2 >${output_folder}/translin.inner
# parallel python country_stat.py --tree ${output_folder}/Sankoff_states.csv.newick --countries $leaf_states_file --max_distance {} \
# --entry_nodes_file ${output_folder}/translin.inner --include_inner --output ${output_folder}/country_stat.include_inner.inner.{} ::: 0 1 2 3

# python entry_parents_withdist.py --tree ${output_folder}/Sankoff_states.csv.newick --steps 1 \
# --entries_file ${translin_file}.entries --output ${output_folder}/translin.parents.1.withdist
# tail -n +2 ${output_folder}/translin.parents.1.withdist | cut -d, -f2  >${output_folder}/translin.parents.1
# parallel python country_stat.py --tree ${output_folder}/Sankoff_states.csv.newick --countries $leaf_states_file --max_distance {} \
# --entry_nodes_file ${output_folder}/translin.parents.1 --include_inner --output ${output_folder}/country_stat.include_inner.parents.1.{} ::: 0 1 2 3

## just test that country_stat works as previously
parallel python country_stat.py --tree ${output_folder}/Sankoff_states.csv.newick --countries $leaf_states_file --max_distance {} \
--entry_nodes_file ${translin_file}.entries --output ${output_folder}/country_stat.test.{} ::: 0 1 2 3

# xvfb-run python -u strip_tree_on_closest.py --tree ${output_folder}/Sankoff_states.csv.newick --melted_foreign ${variants_file}_foreign \
# --dates $leaf_dates_file --entry_nodes_file $variants_file --output ${variants_folder}/stripped_tree --scale 10

## only 3 best entries:
# xvfb-run python -u strip_tree_on_closest.py --tree ${output_folder}/Sankoff_states.csv.newick --melted_foreign ${variants_file}_foreign \
# --dates $leaf_dates_file --entry_nodes_file ${variants_file}_3best --output ${variants_folder}/stripped_tree --scale 10


