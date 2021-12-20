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
mkdir -p $pangolin_output_folder
## when new data is downloaded, cd to previous $pangolin_input_folder and:
## ls *.pangout | xargs -I {} cp {} $pangolin_input_folder/{}
cd $pangolin_input_folder
# conda activate pangolin
## pangolin --update
#find . -regex ".*\.fa[a-z]*" | parallel -j10 pangolin {} --outfile {}.pangout 
## if something goes wrong with pangolin, run find_missing.sh to find corrupted(?) input files ( output to >test_nonexisting) (in /export/home/popova/workspace/covid/data/raw/raw_rus_fasta_29sep)
## cat /export/home/popova/workspace/covid/data/raw/raw_rus_fasta_30sep/test_nonexisting | parallel -j10 pangolin {} --outfile {}.pangout 
#conda deactivate
echo $(pwd)
#grep -n "passed_qc" *.pangout >$rus_pangolined
## something strange happened last time during repangolining, and there are some old pangout files like chuck000001.fa.pangout 
## and also new files like chuck000001.fasta.pangout. Hence instead of just grepping:
find . -name '*.pangout' -mtime 0 -type f -exec grep -n "passed_qc" {} + >$rus_pangolined
sed -i -r "s/.fa[a-z]*.pangout:[0-9]:/,/g" $rus_pangolined
sed -i -r "s!^\.\/!!g" $rus_pangolined
