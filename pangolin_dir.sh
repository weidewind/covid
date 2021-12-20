#!/usr/bin/bash
while getopts c: option
do
		case "${option}" in
			c) config_file=${OPTARG};;
		esac
done

source $config_file
source $conda_config
set -e  

mkdir -p $pangolin_output_folder

# cd $pangolin_input_folder
# conda activate pangolin
# find . -name "*.fasta" | parallel pangolin {} --outfile {}.pangout 
# conda deactivate
# grep -n "passed_qc" *.fasta.pangout >${pangolin_output_folder}/pangout.all
# sed -i "s/.fasta.pangout:2:/,/g" ${pangolin_output_folder}/pangout.all

# cd $script_folder
# python refluff.py --tree ${tree} --duplicates $rus_to_rus_duplicates --output ${tree}.rus_refluffed #add to tree all rus  duplicates of rus strains
# python leaves_without_pangolin_annotation.py --countries $leaf_states_file --tree ${tree}.rus_refluffed --pangolined ${pangolin_output_folder}/pangout.all --output ${pangolin_output_folder}/rus.nopango.idlist

# cd $pangolin_output_folder


# conda activate sarscov2phylo
# echo "faFiltering.."
# faFilter -namePatList=rus.nopango.idlist $unmasked_gisaid_msa rus.nopango.msa
# conda deactivate
# conda activate pangolin
# pangolin rus.nopango.msa --outfile rus.nopango.pangolined
# conda deactivate
# conda activate base

# sed -i -r "s/(EPI_ISL_[0-9]+),/\1,,/g" rus.nopango.pangolined
# sed -i -r "s/taxon,lineage/taxon,stub,lineage/" rus.nopango.pangolined
# ## !! in notepad change (EPI_ISL_[0-9]+), to $1,, nextstrain.pangolined

# cp rus.nopango.pangolined pangolined.all
# cat pangout.all.withdupls >>pangolined.all
# cd $script_folder
python leaves_without_pangolin_annotation.py --countries $leaf_states_file --tree ${tree}.rus_refluffed --pangolined ${pangolin_output_folder}/pangolined.all --output ${pangolin_output_folder}/rus.nopango.test.idlist

##cd ${pangolin_output_folder}
##cat nopango.test.idlist | xargs -I{} grep {} nextstrain.pangolined >>pangolined.all

