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

##seqlength=$( grep ">" -n -m 2 $rus_fasta | tail -1 | cut -d":" -f1 | xargs -I{} head -{} $rus_fasta | head -n -1 | tail -n +2 | tr -d '\n' | wc -c )
## python zero_nodes_debugger.py --tree $tree --output $tree --zero 7.919e-07 # if refluffing was done with a bugged script which adds a new node even if branch length is zero (zero length branches can have length>0, hence --zero option)
# echo "Sequence length is ${seqlength}"
# echo "Will assume that branch lengths are given in ${branchlength}"

# start=`date +%s`
# echo "Metadata conversion.."
# python meta_to_states.py --tree ${tree} --output $leaf_states_file # id	[2 if Russia 1 otherwise]	Country
# python meta_to_dates.py --tree ${tree} --output $leaf_dates_file
# python meta_to_districts.py --tree ${tree} --output $leaf_districts_file
# python meta_to_origins.py --tree ${tree} --output $leaf_origins_file
# end=`date +%s`
# runtime=$( echo "$end - $start" | bc -l )
# echo "metadata converting done in "$((runtime/60))" min "

# echo "State reconstruction.."
# start=`date +%s`
# Rscript Sankoff_asr.R --tree ${tree} --states $leaf_states_file --output $all_states_file
# end=`date +%s`
# runtime=$( echo "$end - $start" | bc -l )
# echo "states reconstructed in "$((runtime/60))" min "

# echo "Searching for transmission lineages.." 
# start=`date +%s`
# python find_transmission_lineages.py --export_rule symm --print_broader_subtree -1 --tree ${all_states_file}.newick --states ${all_states_file}.probs  --countries $leaf_states_file --output ${translin_file} >${translin_file}.log 2>>${translin_file}.errlog
# end=`date +%s`
# runtime=$( echo "$end - $start" | bc -l )
# echo "transmission lineages found in "$((runtime/60))" min "

# echo "Breaking branhces.."
# python entry_branch_breaker.py --tree ${all_states_file}.newick --translin_file $translin_file
 
translin_file=${translin_file}.broken

# echo "Collecting country statistics.."
# start=`date +%s`
# parallel python -u country_stat.py --tree ${all_states_file}.newick.broken --countries $leaf_states_file --entry_nodes_file ${translin_file}.entries \
# --max_distance {} --output ${country_stat_file}.{} >${country_stat_file}.logger.{} ::: 1 2 3 4 5 6
# end=`date +%s`
# runtime=$( echo "$end - $start" | bc -l )
# echo "country stats collected in "$((runtime/60))" min "


# timing ${translin_file}.entries
# timing ${translin_file}.exports


# # stats for dates inside transmission rus_pangolined_from_tree

# echo "Cluster timing.."
# start=`date +%s`
# python cluster_timing.py --transmission_lineages_file $translin_file --dates $leaf_dates_file --output ${output_folder}/cluster_strains.dates_stats
# end=`date +%s`
# runtime=$( echo "$end - $start" | bc -l )
# echo "Cluster timing performed in "$((runtime/60))" min "


# echo "Correcting dates.."

# start=`date +%s`
# for i in $( seq 1 $((threads+1)) )
# do
	# ( Rscript analyze_dates.R  --entry_dates_file ${translin_file}.entries.part${i}.dates_stats --branchlength $branchlength --seqlength $seqlength --output_folder $output_folder ) &
# done
# wait
# cat ${translin_file}.entries.part1.dates_stats.corrected.csv >${translin_file}.entries.dates_stats.corrected.csv
# for i in $( seq 2 $((threads+1)) )
	# do
	# tail -n +2 ${translin_file}.entries.part${i}.dates_stats.corrected.csv >>${translin_file}.entries.dates_stats.corrected.csv
	# done

# echo "Melting entries data.."
# Rscript melt_entries.R  --pangolined_file $rus_pangolined_from_tree  \
 # --entry_strains_file ${translin_file} --leaf_dates_file $leaf_dates_file \
 # --output ${output_folder}/melted_entries.withdates.csv --dates_stats_file ${translin_file}.entries.dates_stats \
 # --entry_dates_file ${translin_file}.entries.dates_stats.corrected.csv

# echo "Drawing plots.."

thr=0.33
# echo ""



# Rscript entry_probability.R --melted_entries_file ${output_folder}/melted_entries.withdates.csv --iterations 1000 --coverage_threshold $thr --output ${output_folder}/entry_freq.coverage_${thr}.csv

# # todo new version:
# Rscript merge_data_for_plots.R --entry_weight_file ${output_folder}/entry_freq.coverage_${thr}.csv \
# --all_pangolined_file $all_pangolined --dates_stats_file ${translin_file}.entries.dates_stats \
# --corrected_dates_file ${translin_file}.entries.dates_stats.corrected.csv \
# --entry_dates_file ${translin_file}.entries.dates_stats --pangolined_file $rus_pangolined_from_tree \
# --cluster_dates_file ${output_folder}/cluster_strains.dates_stats --entry_strains_file ${translin_file} \
# --important_lineages $important_lineages --output ${output_folder}/${tag}_data.csv --coverage $thr

# echo "Drawing all lineages:"
# Rscript rug_plots.R --data_file ${output_folder}/${tag}_data.csv \
# --melted_entries_file ${output_folder}/melted_entries.withdates.csv --districts ${output_folder}/leaf_districts.csv \
# --output_folder $output_folder --tag $tag --coverage $thr --important_lineages $important_lineages

# echo "Drawing selected lineages separately:"
# parallel Rscript rug_plots.R --data_file ${output_folder}/${tag}_data.csv --lineage {} \
# --melted_entries_file ${output_folder}/melted_entries.withdates.csv --districts ${output_folder}/leaf_districts.csv \
# --highlight_secondary \
# --output_folder $output_folder --tag $tag ::: "delta" "AY.122" "B.1.1.397" "B.1.1.397" "B.1.1.317" "BA" "B.1.1.523" "AT.1" 

# parallel Rscript rug_plots.R --data_file ${output_folder}/${tag}_data.csv --lineage {} \
# --melted_entries_file ${output_folder}/melted_entries.withdates.csv --districts ${output_folder}/leaf_districts.csv \
# --output_folder $output_folder --tag $tag ::: "delta" "AY.122" "B.1.1.397" "B.1.1.397" "B.1.1.317" "BA" "B.1.1.523" "AT.1"

# Rscript rug_plots.R --data_file ${output_folder}/${tag}_data.csv --lineage "delta" --sublineage "AY.122" \
# --melted_entries_file ${output_folder}/melted_entries.withdates.csv --districts ${output_folder}/leaf_districts.csv \
# --output_folder $output_folder --tag $tag 


# ###
# # Rscript rug_plots.R --melted_entries_file ${output_folder}/melted_entries.withdates.csv \
# # --data_file ${output_folder}/${tag}_data.csv --districts ${output_folder}/leaf_districts.csv \
# # --output_folder $output_folder --tag ${tag}_melttest --coverage $thr

echo "Pangolin stat plots.." 
Rscript pangolin_stat.R --input $rus_pangolined_from_tree --dates $leaf_dates_file --lineages_list $important_lineages --outfolder $output_folder --title $tag
Rscript pangolin_stat_between.R --input $rus_pangolined_from_tree --dates $leaf_dates_file --lineages_list $important_lineages --outfolder $output_folder --from 2021-07-01 --to 2021-09-15 
Rscript infections_chart.R --input $infections --outfolder $output_folder


echo "Proliferation plots.."
window=14
Rscript proliferation_multilin.R  --window $window --melted_entries_file ${output_folder}/melted_entries.withdates.csv \
--important_lineages $important_lineages --output ${output_folder}/strains_to_entry_deltaggr
Rscript proliferation_multilin.R  --window $window --melted_entries_file ${output_folder}/melted_entries.withdates.csv \
--important_lineages $important_lineages_delta --output ${output_folder}/delta_strains_to_entry --restrict

Rscript clade_size.R --input ${output_folder}/transmission_lineages.withduplicates.outstats --output_folder $output_folder

 end=`date +%s`
 runtime=$( echo "$end - $start" | bc -l )
 echo "Fortnight and rug plots drawn in "$((runtime/60))" min "

