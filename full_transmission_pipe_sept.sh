#!/usr/bin/bash
while getopts c: option
do
		case "${option}" in
			c) config_file=${OPTARG};;
		esac
done

source $config_file

mkdir -p $output_folder
cp $config_file ${output_folder}/$config_file

set -e # exit if any command exits with non-zero status

# #seqlength=$( grep ">" -n -m 2 $rus_fasta | tail -1 | cut -d":" -f1 | xargs -I{} head -{} $rus_fasta | head -n -1 | tail -n +2 | tr -d '\n' | wc -c )
# #python zero_nodes_debugger.py --tree $tree --output $tree --zero 7.919e-07 # if refluffing was done with a bugged script which adds a new node even if branch length is zero (zero length branches can have length>0, hence --zero option)
echo "Sequence length is ${seqlength}"
echo "Will assume that branch lengths are given in ${branchlength}"

start=`date +%s`
echo "Metadata conversion.."
python meta_to_states.py --tree ${tree} --output $leaf_states_file # id	[2 if Russia 1 otherwise]	Country
python meta_to_dates.py --tree ${tree} --output $leaf_dates_file
end=`date +%s`
runtime=$( echo "$end - $start" | bc -l )
echo "metadata converting done in "$((runtime/60))" min "

echo "State reconstruction.."
start=`date +%s`
Rscript Sankoff_asr.R --tree ${tree} --states $leaf_states_file --output $all_states_file
end=`date +%s`
runtime=$( echo "$end - $start" | bc -l )
echo "states reconstructed in "$((runtime/60))" min "

echo "Searching for transmission lineages.." 
start=`date +%s`
python find_transmission_lineages.py --print_broader_subtree -1 --tree ${all_states_file}.newick --states ${all_states_file}.probs  --countries $leaf_states_file --output $translin_file >${translin_file}.log 2>>${translin_file}.errlog
end=`date +%s`
runtime=$( echo "$end - $start" | bc -l )
echo "transmission lineages found in "$((runtime/60))" min "

echo "Collecting country statistics.."
start=`date +%s`
parallel python -u country_stat.py --tree ${all_states_file}.newick --countries $leaf_states_file --entry_nodes_file ${translin_file}.entries \
--max_distance {} --output ${country_stat_file}.{} >${country_stat_file}.logger.{} ::: 1 5 6
end=`date +%s`
runtime=$( echo "$end - $start" | bc -l )
echo "country stats collected in "$((runtime/60))" min "

# Internal nodes dating
# trying to get dates of foreign and russian strains, which are closest to entry node (entry node dating for the poor)
echo "Foreign strains timing.."
start=`date +%s`
#split file with entries into chuncks
entry_num=$(grep -c -v "^$" ${translin_file}.entries)
echo "there are $entry_num lines in ${translin_file}.entries file"

chunk_size=$((entry_num/threads))
filenum=$((threads+1))
echo "chunk size in $chunk_size lines"
for i in $( seq 1 $threads )
	do
	echo "head -$((i*chunk_size)) ${translin_file}.entries | tail -$chunk_size >${translin_file}.entries.part${i}"
	head -$((i*chunk_size)) ${translin_file}.entries | tail -$chunk_size >${translin_file}.entries.part${i}
	done
tail -$((entry_num-threads*chunk_size)) ${translin_file}.entries >${translin_file}.entries.part${filenum}

for i in $( seq 1 $((threads+1)) )
do
	( python -u timing.py --tree ${all_states_file}.newick --dates $leaf_dates_file --countries $leaf_states_file  --entry_nodes_file ${translin_file}.entries.part${i} --output ${translin_file}.part${i} >${translin_file}.part${i}.dates_stats.log ) &
done
wait

cat ${translin_file}.part1.dates_stats >${translin_file}.dates_stats
for i in $( seq 2 $((threads+1)) )
	do
	tail -n +2 ${translin_file}.part${i}.dates_stats >>${translin_file}.dates_stats
	#unlink ${translin_file}.part${i}.dates_stats
	#unlink ${translin_file}.part${i}.dates_stats.log
	done

end=`date +%s`
runtime=$( echo "$end - $start" | bc -l )
echo "Entry timing performed in "$((runtime/60))" min "

# stats for dates inside transmission rus_withdupls_pangolined

echo "Cluster timing.."
start=`date +%s`
python cluster_timing.py --transmission_lineages_file $translin_file --dates $leaf_dates_file --output ${output_folder}/cluster_strains.dates_stats
end=`date +%s`
runtime=$( echo "$end - $start" | bc -l )
echo "Cluster timing performed in "$((runtime/60))" min "


echo "Correcting dates.."

start=`date +%s`
for i in $( seq 1 $((threads+1)) )
do
	( Rscript analyze_dates.R  --entry_dates_file ${translin_file}.part${i}.dates_stats --branchlength $branchlength --seqlength $seqlength --output_folder $output_folder ) &
done
wait
cat ${translin_file}.part1.dates_stats.corrected.csv >${translin_file}.dates_stats.corrected.csv
for i in $( seq 2 $((threads+1)) )
	do
	tail -n +2 ${translin_file}.part${i}.dates_stats.corrected.csv >>${translin_file}.dates_stats.corrected.csv
	done

echo "Melting entries data.."
Rscript melt_entries.R  --pangolined_file $rus_withdupls_pangolined  \
 --entry_strains_file ${translin_file} --leaf_dates_file $leaf_dates_file \
 --output ${output_folder}/melted_entries.withdates.csv --dates_stats_file ${translin_file}.dates_stats \
 --entry_dates_file ${translin_file}.dates_stats.corrected.csv

echo "Drawing plots.."

thr=0.33
Rscript entry_probability.R --melted_entries_file ${output_folder}/melted_entries.withdates.csv --iterations 1000 --coverage_threshold $thr --output ${output_folder}/entry_freq.coverage_${thr}.csv

## todo new version:
Rscript merge_data_for_plots.R --entry_weight_file ${output_folder}/entry_freq.coverage_${thr}.csv \
--corrected_dates_file ${translin_file}.dates_stats.corrected.csv \
--entry_dates_file ${translin_file}.dates_stats --pangolined_file $rus_withdupls_pangolined \
--cluster_dates_file ${output_folder}/cluster_strains.dates_stats --entry_strains_file ${translin_file} \
--important_lineages $important_lineages --output ${output_folder}/${tag}_data.csv --coverage $thr

Rscript rug_plots.R --data_file ${output_folder}/${tag}_data.csv \
--melted_entries_file ${output_folder}/melted_entries.withdates.csv --districts ${output_folder}/leaf_districts.csv \
--output_folder $output_folder --tag $tag --coverage $thr

Rscript rug_plots.R --data_file ${output_folder}/${tag}_data.csv --lineage "delta" \
--melted_entries_file ${output_folder}/melted_entries.withdates.csv --districts ${output_folder}/leaf_districts.csv \
--output_folder $output_folder --tag $tag 

Rscript rug_plots.R --data_file ${output_folder}/${tag}_data.csv --lineage "B.1.1.397" \
--melted_entries_file ${output_folder}/melted_entries.withdates.csv --districts ${output_folder}/leaf_districts.csv \
--output_folder $output_folder --tag $tag 

Rscript rug_plots.R --data_file ${output_folder}/${tag}_data.csv --lineage "B.1.1.317" \
--melted_entries_file ${output_folder}/melted_entries.withdates.csv --districts ${output_folder}/leaf_districts.csv \
--output_folder $output_folder --tag $tag 

# Rscript rug_plots.R --data_file ${output_folder}/${tag}_data.csv --lineage "AT.1" \
# --melted_entries_file ${output_folder}/melted_entries.withdates.csv --districts ${output_folder}/leaf_districts.csv \
# --output_folder $output_folder --tag $tag 
##

Rscript rug_plots.R --melted_entries_file ${output_folder}/melted_entries.withdates.csv \
--entry_weight_file ${output_folder}/entry_freq.coverage_${thr}.csv --corrected_dates_file ${translin_file}.dates_stats.corrected.csv \
--entry_dates_file ${translin_file}.dates_stats --pangolined_file $rus_withdupls_pangolined \
--cluster_dates_file ${output_folder}/cluster_strains.dates_stats --entry_strains_file ${translin_file} \
--important_lineages $important_lineages --output_folder $output_folder --tag ${tag}_melttest --coverage $thr

echo "Pangolin stat plots.." 
Rscript pangolin_stat.R --input $rus_withdupls_pangolined --dates $leaf_dates_file --lineages_list $important_lineages --outfolder $output_folder --title $tag
Rscript pangolin_stat_between.R --input $rus_withdupls_pangolined --dates $leaf_dates_file --lineages_list $important_lineages --outfolder $output_folder --from 2021-07-01 --to 2021-09-15 
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

