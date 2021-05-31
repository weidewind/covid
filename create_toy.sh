#!/usr/bin/bash

cd ~/workspace/covid/data

refid="EPI_ISL_406798"
faOneRecord raw/msa_2021-03-22/2021-03-22_unmasked.fa $refid >toy/reference.fa

head -24500  raw/msa_2021-03-22/2021-03-22_unmasked.fa >toy/global.msa # first 50 seq
cat toy/reference.fa >>toy/global.msa

tail -n +1107644 russian/alignment/all_rus_nodupl_mafft.aln >toy/rus.aln # last 10 seq (because why not)

cd toy

iqtree -s global.msa -m GTR+I+G

#toy test add_strains_to_tree.sh
#./add_strains_to_tree.sh -t  ~/workspace/covid/data/toy/global.msa.treefile -f ~/workspace/covid/data/toy/rus.aln -g ~/workspace/covid/data/toy/global.msa -n 3


