#!/usr/bin/bash

DIR="/export/home/popova/workspace/covid/output/transmission_lineages/gennady/far20"
for t in ${DIR}/transmission_lineages.withduplicates.out_entry_*_5_steps_up_newick; do
	treetime --keep-root --sequence-length 29903 --tree $t --dates ${DIR}/leaf_dates.csv.cleaned.min  --outdir ${DIR}/separate >>${DIR}/separate/treetime.logs
done