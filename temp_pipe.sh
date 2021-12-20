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
echo "Sequence length is ${seqlength}"
echo "Will assume that branch lengths are given in ${branchlength}"

echo "Collecting country statistics.."
start=`date +%s`
parallel python -u country_stat.py --tree ${all_states_file}.newick --states ${all_states_file}.probs --countries $leaf_states_file --entry_nodes_file ${translin_file}.entries --max_distance {} --output ${country_stat_file}.{} >${country_stat_file}.logger.{} ::: 0.00025 0.00021
end=`date +%s`
runtime=$( echo "$end - $start" | bc -l )
echo "country stats collected in "$((runtime/60))" min "