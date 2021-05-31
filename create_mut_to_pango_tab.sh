#!/usr/bin/bash

cd ../data/raw
grep "Russia" msa_0319/msa_0319.fasta | sed 's/^>hCoV-19\///' |cut -d "|" -f '1 2' >id_to_epiisl.msa_0319.txt
# немножко экселя в russian/meta_all.xls

cd ../../scripts
conda activate sarscov2phylo
faFilter -v -namePatList=data/raw/our_ids_at_gisaid_0319.txt data/russian/alignment/all_rus_mafft.aln data/russian/alignment/all_rus_nodupl_mafft.aln
conda deactivate

conda activate pangolin
pangolin data/russian/alignment/all_rus_nodupl_mafft.aln --outfile data/russian/alignment/all_rus_nodupl.pangolin

python scripts/mutation_caller.py --fasta data/russian/alignment/all_rus_nodupl_mafft.aln --output data/russian/alignment/all_rus.mutations --mutations data/russian/alignment/rpn.risky
python scripts/pango_pandas.py --pango_file data/russian/alignment/all_rus_nodupl.pangolin --mut_file data/russian/alignment/all_rus.mutations.table.csv --output data/russian/alignment/all_rus_merged.csv

# lets do the same with a nice protein alignment of spike protein
# cut and realign, since inserted gaps are untranslatable now (ATTENTION magic constants inside)
python scripts/cut_alignment.py --fasta data/russian/alignment/all_rus_nodupl_mafft.aln --output data/russian/alignment/all_rus_nodupl_spike.aln --start 21563 --end 25384
#todo conventional fasta without gaps -> translate (catch and exclude non-translatable) -> align protein


