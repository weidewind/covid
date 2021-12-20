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



## Pangoline all russian strains present in $tree AND in $unmasked_gisaid_msa, but not in $previously_pangolined. 
## $previously_pangolined must have a header starting from "taxon,stub,lineage"

# python meta_to_states.py --tree $tree --output $leaf_states_file
# python leaves_without_pangolin_annotation.py --tree $tree --countries $leaf_states_file --pangolined $previously_pangolined --output ${pangolin_output_folder}/rus.nopango.idlist
# cd $pangolin_output_folder
# conda activate sarscov2phylo
# echo "faFiltering gisaid strains without pangolin annotation.."
# faFilter -namePatList=rus.nopango.idlist $unmasked_gisaid_msa rus.nopango.msa
# conda deactivate
# conda activate pangolin
# pangolin rus.nopango.msa --outfile rus.nopango.pangolined
# sed -i -r "s/(EPI_ISL_[0-9]+),/\1,,/g" rus.nopango.pangolined
# sed -i -r "s/taxon,lineage/taxon,stub,lineage/" rus.nopango.pangolined
cp rus.nopango.pangolined $new_pangolined
tail -n +2 $previously_pangolined >>$new_pangolined

cd $script_folder
echo $leaf_states_file
echo $new_pangolined
echo $rus_withdupls_pangolined
python select_rus_pangolined.py --tree $tree --states $leaf_states_file --pangolined $new_pangolined --output $rus_withdupls_pangolined

