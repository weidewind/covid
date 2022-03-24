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

## pangolin all our strains (raw data), no filtering applied

mkdir -p $pangolin_output_folder # (export/.../data/munched/pangolin_v..)
## when new data is downloaded and pangolin version is unchanged, copy previous results, if they exist:
if [ -n "$previous_pangolin_output_folder" ]
then
	echo "Copying previously pangolined files from $previous_pangolin_output_folder .."
	cd $previous_pangolin_output_folder
	ls *.pangout | xargs -I {} cp {} ${pangolin_output_folder}/{}
fi

echo "Making a list of fasta files with strains present in the tree, which have no corresponding pangouts.. will print them into ${pangolin_output_folder}/missing.list and copy to ${output_folder}/missing.list"
${script_folder}/find_missing.sh -i $pangolin_input_folder -o $pangolin_output_folder -t $tree >${pangolin_output_folder}/missing.list
cp ${pangolin_output_folder}/missing.list ${output_folder}/missing.list
echo "If non-ascii chars are present in any of these files, change them to N.."
cd $pangolin_input_folder
cat ${pangolin_output_folder}/missing.list | xargs -I {} perl -i -pa -ne 's/[[:^ascii:]]/N/g' {}


conda activate pangolin
## pangolin --update
echo "Pangolining.."
cat ${pangolin_output_folder}/missing.list | parallel -j10 pangolin {} --outfile $pangolin_output_folder/{}.pangout 
## if something goes wrong with pangolin, run find_missing.sh to find corrupted(?) input files ( output to >test_nonexisting) (in /export/home/popova/workspace/covid/data/raw/raw_rus_fasta_29sep)
echo "Checking if some pangouts are still missing.. will write them to ${pangolin_output_folder}/missing.list and try again"
${script_folder}/find_missing.sh -i $pangolin_input_folder -o $pangolin_output_folder -t $tree >${pangolin_output_folder}/missing.errlist
cp ${pangolin_output_folder}/missing.errlist ${output_folder}/missing.errlist
cat ${pangolin_output_folder}/missing.errlist | parallel -j10 pangolin {} --outfile ${pangolin_output_folder}/{}.pangout 
conda deactivate

cd $pangolin_output_folder
echo "Collecting pangout files from $pangolin_output_folder into file $our_pangolined .."
## grep -n "passed_qc" *.pangout >$our_pangolined # Does not work, because Argument list too long
## something strange happened last time during repangolining, and there are some old pangout files like chuck000001.fa.pangout 
## and also new files like chuck000001.fasta.pangout. Hence instead of just grepping:
## find . -name '*.pangout' -mtime 0 -type f -exec grep -n "passed_qc" {} + >$our_pangolined
find . -name '*.pangout' -type f -exec grep -n "passed_qc" {} + >$our_pangolined
sed -i -r "s/.fa[a-z]*.pangout:[0-9]:/,/g" $our_pangolined
sed -i -r "s!^\.\/!!g" $our_pangolined
echo "All done!"
