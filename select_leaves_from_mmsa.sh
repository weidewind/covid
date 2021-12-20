#!/usr/bin/bash

tag="3oct_opt_min_newpango"
conda_config="/export/home/popova/miniconda3/etc/profile.d/conda.sh"
#masked_gisaid_msa=/export/home/popova/workspace/covid/data/raw/mmsa_2021-11-17/2021-11-17_masked.fa
deprecated_masked_gisaid_msa=/export/home/popova/workspace/covid/data/raw/mmsa_2021-09-23/mmsa_2021-09-23/2021-09-23_masked.fa
variants_folder=/export/home/popova/workspace/covid/output/variant_analysis/${tag}
rus_msa=/export/home/popova/workspace/covid/data/raw/raw_rus_fasta_30sep/2021-09-30_all_rus_for_usher_nextstrain.fasta
mao=/export/home/popova/workspace/ksuAncestorPipeline/ancestral_seqs_ML_nucleotide.mao
duplicates=/export/home/popova/workspace/covid/data/munched/usher/3oct_opt/2021-09-30_rus_to_rus_dups_not_too_far.tsv

tree=${variants_folder}/3clades.nozero.bif.uniroot.nosingle.nwk
#mytree=/export/home/popova/workspace/ksuAncestorPipeline/h1pandemic/tree.rooted.nwk
#testtree=${variants_folder}/test.nwk

source $conda_config

# python get_leaves.py --tree $tree --output ${tree}.idlist
# conda activate sarscov2phylo
# faFilter -namePatList=${tree}.idlist $deprecated_masked_gisaid_msa ${tree}.gisaid.fasta
# >${tree}.notfound.idlist
# while IFS="" read -r id || [ -n "$id" ] # otherwise it will skip the trailing line
# do
  # grep "${id}," ${tree}.gisaid.fasta || echo $id >>${tree}.notfound.idlist
# done <${tree}.idlist
# faFilter -namePatList=${tree}.notfound.idlist $rus_msa ${tree}.next.fasta
# cp ${tree}.gisaid.fasta ${tree}.all.fasta
# cat ${tree}.next.fasta >>${tree}.all.fasta
# conda deactivate

## cut tails to get the same length for all sequences - 29880 (what are these polyA at the end of rus msa seqs?)
# python cut_tails.py --fasta ${tree}.all.fasta --length 29880 --output ${tree}.all.cut.fasta


python sub_msa.py --fasta ${tree}.all.cut.fasta --tree $tree --duplicates $duplicates --sites 1048,27527 --output ${tree}.all.2sites.fasta
xvfb-run megacc -a $mao -d ${tree}.all.2sites.fasta -f Fasta -t $tree -o ${tree}.2sites.ancestors

