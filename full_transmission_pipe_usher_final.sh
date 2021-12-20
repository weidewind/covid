#!/usr/bin/bash
while getopts c: option
do
		case "${option}" in
			c) config_file=${OPTARG};;
		esac
done

source $config_file

mkdir -p $output_folder

set -e # exit if any command exits with non-zero status

seqlength=$( grep ">" -n -m 2 $rus_fasta | tail -1 | cut -d":" -f1 | xargs -I{} head -{} $rus_fasta | head -n -1 | tail -n +2 | tr -d '\n' | wc -c )
echo "Sequence length is ${seqlength}"
echo "Branch lengths are given in ${branchlength}"


echo "Refluffing and renaming.."
start=`date +%s`
python refluff.py --tree $rawtree --duplicates $rus_to_rus_duplicates --output ${rawtree}.rus_refluffed #add to tree all rus  duplicates of rus strains
python rename_leaves.py --tree ${rawtree}.rus_refluffed --meta $merged_meta --long_ids $long_ids --output ${tree}.rus_refluffed
echo "Metadata conversion.."
python meta_to_states.py --tree ${tree}.rus_refluffed --output $leaf_states_file # id	[2 if Russia 1 otherwise]	Country
python meta_to_dates.py --tree ${tree}.rus_refluffed --output $leaf_dates_file
end=`date +%s`
runtime=$( echo "$end - $start" | bc -l )
echo "metadata converting and tree refluffing done in "$((runtime/60))" min "


echo "State reconstruction.."
start=`date +%s`
Rscript Sankoff_asr.R --tree ${tree}.rus_refluffed --states $leaf_states_file --output $all_states_file
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
parallel python -u country_stat.py --tree ${all_states_file}.newick --states ${all_states_file}.probs --countries $leaf_states_file --entry_nodes_file ${translin_file}.entries --max_distance {} --output ${country_stat_file}.{} >${country_stat_file}.logger.{} ::: 0 1 2 3 4 5 6
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

# stats for dates inside transmission clusters
echo "Cluster timing.."
start=`date +%s`
python cluster_timing.py --transmission_lineages_file $translin_file --dates $leaf_dates_file --output ${output_folder}/cluster_strains.dates_stats
end=`date +%s`
runtime=$( echo "$end - $start" | bc -l )
echo "Cluster timing performed in "$((runtime/60))" min "

echo "Pangolining russian strains.."
start=`date +%s`
mkdir -p $pangolin_output_folder
source $conda_config
# cd $pangolin_input_folder
# conda activate pangolin
# find . -name "*.fasta" | parallel pangolin {} --outfile {}.pangout 
# conda deactivate
# grep -n "passed_qc" *.fasta.pangout >${pangolin_output_folder}/pangout.all
# sed -i "s/.fasta.pangout:2:/,/g" ${pangolin_output_folder}/pangout.all
cd $script_folder
python leaves_without_pangolin_annotation.py --countries $leaf_states_file --tree ${tree}.rus_refluffed --pangolined ${pangolin_output_folder}/pangout.all --output ${pangolin_output_folder}/rus.nopango.idlist
cd $pangolin_output_folder
conda activate sarscov2phylo
echo "faFiltering rus gisaid strains without pangolin annotation.."
faFilter -namePatList=rus.nopango.idlist $unmasked_gisaid_msa rus.nopango.msa
conda deactivate
conda activate pangolin
pangolin rus.nopango.msa --outfile rus.nopango.pangolined
sed -i -r "s/(EPI_ISL_[0-9]+),/\1,,/g" rus.nopango.pangolined
sed -i -r "s/taxon,lineage/taxon,stub,lineage/" rus.nopango.pangolined
cp rus.nopango.pangolined $rus_withdups_pangolined
cat pangout.all.withdupls >>$rus_withdups_pangolined
cd $script_folder
python leaves_without_pangolin_annotation.py --countries $leaf_states_file --tree ${tree}.rus_refluffed --pangolined ${pangolin_output_folder}/panglined.all --output ${pangolin_output_folder}/rus.nopango.test.idlist
# check that ${pangolin_output_folder}/nopango.test.idlist is empty!
conda deactivate
conda activate base
end=`date +%s`
runtime=$( echo "$end - $start" | bc -l )
echo "Pangolining performed in "$((runtime/60))" min "

echo "Drawing pangolin plots.." 
start=`date +%s`
Rscript pangostat.R --pangolined_file $rus_withdups_pangolined --entry_file ${translin_file} --important_lineages $important_lineages --output_folder $pangolin_output_folder --tag $gentag
end=`date +%s`
runtime=$( echo "$end - $start" | bc -l )
echo "Pangolin plots drawn in "$((runtime/60))" min "


echo "Drawing fortnight and rug plots.."
#instead of #Rscript analyze_dates.R --pangolined_file $rus_withdups_pangolined  --cluster_dates_file ${output_folder}/cluster_strains.dates_stats --entry_strains_file ${translin_file} --entry_dates_file ${translin_file}.dates_stats --important_lineages $important_lineages --output_folder $pangolin_output_folder --tag $gentag
start=`date +%s`
for i in $( seq 1 $((threads+1)) )
do
	( Rscript analyze_dates.R  --entry_dates_file ${translin_file}.part${i}.dates_stats --branchlength $branchlength --seqlength $seqlength --output_folder $pangolin_output_folder --tag $gentag ) &
done
wait
cat ${translin_file}.part1.dates_stats.corrected.${gentag}.csv >${translin_file}.dates_stats.corrected.${gentag}.csv
for i in $( seq 2 $((threads+1)) )
	do
	tail -n +2 ${translin_file}.part${i}.dates_stats.corrected.${gentag}.csv >>${translin_file}.dates_stats.corrected.${gentag}.csv
	done

Rscript rug_plots.R --corrected_dates_file ${translin_file}.dates_stats.corrected.${gentag}.csv --entry_dates_file ${translin_file}.dates_stats --pangolined_file $rus_withdups_pangolined --cluster_dates_file ${output_folder}/cluster_strains.dates_stats --entry_strains_file ${translin_file} --important_lineages $important_lineages --output_folder $pangolin_output_folder --tag $gentag
end=`date +%s`
runtime=$( echo "$end - $start" | bc -l )
echo "Fortnight and rug plots drawn in "$((runtime/60))" min "


mkdir -p $output_folder/trees
start=`date +%s`
xvfb-run python draw_stripped_pangolineages_of_interest.py --lineages_list $important_lineages \
--pangostats_file ${pangolin_output_folder}/${gentag}_pangolin_counts_for_entry.csv --tree ${all_states_file}.newick \
--steps 3 --countries $leaf_states_file --dates $leaf_dates_file \
--pangolined $rus_withdups_pangolined --states ${all_states_file}.probs --entries_file $translin_file \
--output ${output_folder}/trees/${gentag}_entries_of_interest
end=`date +%s`
runtime=$( echo "$end - $start" | bc -l )
echo "stripped pangolineages of interest drawn in "$((runtime/60))" min "


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
