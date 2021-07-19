#!/usr/bin/bash

tag="nextstrain_and_17june_3"
mode="c" #c circular or r rectangular
outtree="nexttree_17june_3_${mode}_large"
rawdatafolder=/export/home/popova/workspace/covid/data/raw/nextstrain
datafolder=/export/home/popova/workspace/covid/data/munched/usher/${tag}
outfolder=/export/home/popova/workspace/covid/output/nextstrain/${tag}
dates_file=${datafolder}/leaf_dates.csv
rus_to_rus_duplicates=${datafolder}/2021-07-01_rus_to_rus_dups.tsv

tree=${datafolder}/uncondensed-final-tree.nh
refluffed_tree=${tree}.refluffed
rus_aln=${datafolder}/2021-07-01_all_rus_for_usher_nextstrain.fasta # 17_06_2021_rus_masked_with_ref_trimmed.aln
nextstrain_aln=${datafolder}/nextstrain_ncov_global_timetree_21_06_202_gisaid_names_ref_first.aln
rus_pangolined=${outfolder}/rus.pangolined
nextstrain_pangolined=${outfolder}/nextstrain.pangolined
all_pangolined=${outfolder}/all.pangolined
all_withdups_pangolined=${outfolder}/pangolined.all

nextmeta=${rawdatafolder}/nextstrain_ncov_global_metadata_21_06_2021.tsv
rusmeta=/export/home/popova/workspace/covid/data/raw/17_06_2021_meta.csv

conda_config="/export/home/popova/miniconda3/etc/profile.d/conda.sh"

mkdir -p $outfolder

# tar -xvf /export/data/popova/covid/raw/mmsa_2021-06-21.tar.xz
# conda activate sarscov2phylo
# cut -f12 nextstrain_ncov_global_metadata_21_06_2021.tsv >idlist
# nohup faFilter -namePatList=idlist mmsa_2021-06-21/2021-06-21_masked.fa nextstrain_align.fa &
# grep -c ">" nextstrain_align.fa
# wc -l idlist

#sed удалить все " из метаданных некстстрейна
# убрать из русского выравнивания то, что уже есть в некстстрейне
# добивать в начало обоих выравниваний правильный референс 



# source $conda_config
# conda deactivate
# conda activate pangolin
# pangolin $rus_aln --outfile $rus_pangolined
# pangolin $nextstrain_aln --outfile $nextstrain_pangolined
# cp $rus_pangolined $all_pangolined
# cat $nextstrain_pangolined >>$all_pangolined
# conda deactivate

#python refluff.py --tree $tree --duplicates $rus_to_rus_duplicates --output $refluffed_tree
#python pango_add_duplicates_anarchy.py --pango_file $all_pangolined --duplicates_file $rus_to_rus_duplicates --output $all_withdups_pangolined

# echo "xvfb-run python draw_nextstrain_tree.py --tree $refluffed_tree --rusmeta $rusmeta --nextmeta $nextmeta --pangolined $all_pangolined --output ${outfolder}/${outtree}"
xvfb-run python draw_nextstrain_tree.py --tree $refluffed_tree --rusmeta $rusmeta --nextmeta $nextmeta --dates $dates_file --pangolined $all_withdups_pangolined --mode $mode --output ${outfolder}/${outtree}
echo "Editing graphics.."
newr=10
if [ "$mode" = "c" ]; then
	newr=45
fi
echo $newr
sed -i 's/<g/\n<g/g' ${outfolder}/${outtree}.svg
sed -i "s/r=\"2\"/r=\"$newr\"/g" ${outfolder}/${outtree}.svg # 
grep -v "circle" ${outfolder}/${outtree}.svg | head -n -1 >${outfolder}/${outtree}_edited.svg
grep "circle" ${outfolder}/${outtree}.svg >>${outfolder}/${outtree}_edited.svg
tail -1 ${outfolder}/${outtree}.svg >>${outfolder}/${outtree}_edited.svg
#perl -pi -e 's/viewBox(.{2})([0-9\.]+) ([0-9\.]+) ([0-9\.]+) ([0-9\.]+)/viewBox.$1.($2-800)." ".($3-800)." ".($4+1600)." ".($5+1600)/ge' ${outfolder}/${outtree}_edited.svg



python convert_svg_to_pdf.py --input ${outfolder}/${outtree}_edited.svg --output ${outfolder}/${outtree}_edited.pdf

