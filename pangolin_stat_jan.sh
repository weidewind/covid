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

## Pangolin gisaid strains present in the tree which has no pangolin annotation in our.pangolined
cd $script_folder
echo "Looking for strains present in the tree which has no pangolin annotation in our.pangolined.."
python leaves_without_pangolin_annotation.py --tree $tree --pangolined ${our_pangolined},${nopango_pangolined} --output ${output_folder}/nopango.idlist
cd $output_folder
conda activate sarscov2phylo
echo "faFiltering gisaid strains without pangolin annotation.."
faFilter -namePatList=nopango.idlist $unmasked_gisaid_msa nopango.msa
conda deactivate
conda activate pangolin
echo "pangolining.."
pangolin nopango.msa --outfile ${nopango_pangolined}.temp
conda deactivate
conda activate base
echo "editing and merging.."
sed -i -r "s/(EPI_ISL_[0-9]+),/\1,,/g" ${nopango_pangolined}.temp
sed -i -r "s/taxon,lineage/taxon,stub,lineage/" ${nopango_pangolined}.temp
if [ -f $nopango_pangolined ]
then
	tail -n +2 ${nopango_pangolined}.temp >>$nopango_pangolined
else
	mv ${nopango_pangolined}.temp $nopango_pangolined
fi
cat $nopango_pangolined >$all_pangolined
cat $our_pangolined >>$all_pangolined
cd $script_folder

## check if some strains still has no pangolin annotation
echo "checking if any strains are still unpangolined.. will output it to "
echo ${output_folder}/nopango.idlist.errlist
python leaves_without_pangolin_annotation.py --tree $tree --pangolined $all_pangolined --output ${output_folder}/nopango.idlist.errlist

## select russian strains present in the tree and draw plots
# python meta_to_states.py --tree $tree --output $leaf_states_file
echo "drawing plots.."
python select_rus_pangolined.py --tree $tree --states $leaf_states_file --pangolined $all_pangolined --output $rus_pangolined_from_tree
Rscript pangolin_stat.R --input $rus_pangolined_from_tree --dates $leaf_dates_file --lineages_list $important_lineages --outfolder $output_folder --title $tag