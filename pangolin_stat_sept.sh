#!/usr/bin/bash
while getopts c: option
do
		case "${option}" in
			c) config_file=${OPTARG};;
		esac
done

source $config_file
source $conda_config
conda deactivate
conda activate base

# python meta_to_states.py --tree $tree --output $leaf_states_file
# python select_rus_pangolined.py --tree $tree --output $rus_pangolined_from_tree --states $leaf_states_file --pangolined $all_pangolined
# head -1 $all_pangolined >${outfolder}/pango.temp
# cat $rus_pangolined_from_tree >>${outfolder}/pango.temp
# mv ${outfolder}/pango.temp $rus_pangolined_from_tree
Rscript pangolin_stat.R --input $rus_pangolined_from_tree --dates $leaf_dates_file --lineages_list $important_lineages --outfolder $outfolder --title $tag