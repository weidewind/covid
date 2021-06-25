#!/usr/bin/bash

tag="may31"
output_folder="../data/munched/gennady/${tag}"
rus_unique_fasta="${output_folder}/rus_unique.fasta"
rus_to_rus_duplicates="${output_folder}/21_05_2021_rus_to_rus_dups.tsv"
unique_pangolined="${output_folder}/rus_unique.pangolined"
pangolined_with_duplicates="${output_folder}/rus.pangolined.withduplicates"
lineages_of_interest="/export/home/popova/workspace/covid/data/russian/lineages_of_interest.list"

conda activate pangolin
pangolin $rus_unique_fasta --outfile $unique_pangolined
conda deactivate
python --pango_file $unique_pangolined --duplicates_file $rus_to_rus_duplicates --output $pangolined_with_duplicates

#pangolin_stat.Rmd #todo make it a script

