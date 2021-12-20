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

# cd $script_folder
# python leaves_without_pangolin_annotation.py --countries $leaf_states_file --tree $tree --pangolined $rus_pangolined --output ${pangolin_output_folder}/rus.nopango.idlist
# cd $pangolin_output_folder
# conda activate sarscov2phylo
# echo "faFiltering gisaid strains without pangolin annotation.."
# faFilter -namePatList=rus.nopango.idlist $unmasked_gisaid_msa rus.nopango.msa
# conda deactivate
# conda activate pangolin
# pangolin rus.nopango.msa --outfile rus.nopango.pangolined
# sed -i -r "s/(EPI_ISL_[0-9]+),/\1,,/g" rus.nopango.pangolined
# sed -i -r "s/taxon,lineage/taxon,stub,lineage/" rus.nopango.pangolined
# cp rus.nopango.pangolined $all_rus_pangolined
# cat $rus_pangolined >>$all_rus_pangolined

cd $script_folder
python select_rus_pangolined.py --tree $tree --output $rus_pangolined_from_tree_shortcut --states $leaf_states_file --pangolined $all_rus_pangolined
Rscript pangolin_stat.R --input $rus_pangolined_from_tree_shortcut --dates $leaf_dates_file --lineages_list $important_lineages --outfolder $output_folder --title $tag