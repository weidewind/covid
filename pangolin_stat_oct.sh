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

# We already have $rus_pangolined from pangolin_dir_sept.sh. Now:
# TODO!
#1) select from the tree russian sequences without annotation, pangolin them and merge with $rus_pangolined to get $all_pangolined
#2) select russian sequences present in the tree to cut them from $all_pangolined to get $rus_withdups_pangolined

cd $script_folder
#--countries $leaf_states_file
python leaves_without_pangolin_annotation.py --tree $tree --pangolined $rus_pangolined --output ${pangolin_output_folder}/nopango.idlist
cd $pangolin_output_folder
conda activate sarscov2phylo
echo "faFiltering gisaid strains without pangolin annotation.."
faFilter -namePatList=nopango.idlist $unmasked_gisaid_msa nopango.msa
conda deactivate
conda activate pangolin
pangolin nopango.msa --outfile nopango.pangolined
sed -i -r "s/(EPI_ISL_[0-9]+),/\1,,/g" nopango.pangolined
sed -i -r "s/taxon,lineage/taxon,stub,lineage/" nopango.pangolined
cp nopango.pangolined $all_pangolined
cat $rus_pangolined >>$all_pangolined
cd $script_folder

python select_rus_pangolined.py --tree $tree --output $rus_pangolined_from_tree --states $leaf_states_file --pangolined $all_pangolined
Rscript pangolin_stat.R --input $rus_pangolined_from_tree --dates $leaf_dates_file --lineages_list $important_lineages --outfolder $output_folder --title $tag