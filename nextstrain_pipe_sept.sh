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


source $conda_config
mkdir -p outfolder
cp $config_file ${outfolder}/$config_file

# # tar -xvf /export/data/popova/covid/raw/mmsa_2021-06-21.tar.xz
# # conda activate sarscov2phylo
# # cut -f12 nextstrain_ncov_global_metadata_21_06_2021.tsv >idlist
# # nohup faFilter -namePatList=idlist mmsa_2021-06-21/2021-06-21_masked.fa nextstrain_align.fa &
# # grep -c ">" nextstrain_align.fa
# # wc -l idlist

# # sed удалить все " из метаданных некстстрейна
# # убрать из русского выравнивания то, что уже есть в некстстрейне
# # добивать в начало обоих выравниваний правильный референс 

# conda deactivate
# conda activate pangolin
# pangolin $nextstrain_aln --outfile $nextstrain_pangolined
# conda deactivate
# sed -r "s/(EPI_ISL_[0-9]+),/\1,,/g" $nextstrain_pangolined >${nextstrain_pangolined}.edited
# sed -i "s/taxon,lineage/taxon,stub,lineage/g" ${nextstrain_pangolined}.edited

# cat ${nextstrain_pangolined}.edited >$all_pangolined
# cat $rus_pangolined >>$all_pangolined

# conda activate base
# python leaves_without_pangolin_annotation.py --tree $tree --pangolined $all_pangolined --output ${outfolder}/nopango.idlist
# conda activate sarscov2phylo
# echo "faFiltering.."
# faFilter -namePatList=${outfolder}/nopango.idlist $unmasked_gisaid_msa ${outfolder}/nopango.msa
# conda deactivate
# conda activate pangolin
# pangolin ${outfolder}/nopango.msa --outfile ${outfolder}/nopango.pangolined
# conda deactivate
# conda activate base
# sed -r "s/(EPI_ISL_[0-9]+),/\1,,/g" ${outfolder}/nopango.pangolined >${outfolder}/nopango.pangolined.edited
# cat ${outfolder}/nopango.pangolined.edited | tail -n +2 >>$all_pangolined

#no nextstrain metadata is available at nextstrain.org now, so we use gisaid meta instead
# python meta_to_dates.py --tree $tree --output $leaf_dates_file
outtree=${outtree}_unified_colors_test_overlap_aligned
xvfb-run python draw_nextstrain_tree_sept.py --tree $tree  --states $leaf_states_file --dates $leaf_dates_file --important_lineages $important_lineages --pangolined $all_pangolined --mode $mode --output ${outfolder}/${outtree} --print_delta
echo "Editing graphics.."
newr=10
if [ "$mode" = "c" ]; then
	newr=70
fi
echo $newr
sed 's/<g/\n<g/g' ${outfolder}/${outtree}.svg >${outfolder}/${outtree}_temp.svg
sed -i "s/r=\"2\"/r=\"$newr\"/g" ${outfolder}/${outtree}_temp.svg
grep -v "circle" ${outfolder}/${outtree}_temp.svg | head -n -1 >${outfolder}/${outtree}_edited.svg
grep "circle" ${outfolder}/${outtree}_temp.svg >>${outfolder}/${outtree}_edited.svg
tail -1 ${outfolder}/${outtree}_temp.svg >>${outfolder}/${outtree}_edited.svg
perl -pi -e 's/viewBox(.{2})([0-9\.]+) ([0-9\.]+) ([0-9\.]+) ([0-9\.]+)/viewBox.$1.($2-800)." ".($3-800)." ".($4+1600)." ".($5+1600)/ge' ${outfolder}/${outtree}_edited.svg

python convert_svg_to_pdf.py --input ${outfolder}/${outtree}_edited.svg --output ${outfolder}/${outtree}_edited.pdf

newr=20
sed 's/<g/\n<g/g'  ${outfolder}/${outtree}_delta.svg >${outfolder}/${outtree}_delta_temp.svg
sed -i "s/r=\"2\"/r=\"$newr\"/g" ${outfolder}/${outtree}_delta_temp.svg # 
grep -v "circle" ${outfolder}/${outtree}_delta_temp.svg | head -n -1 >${outfolder}/${outtree}_delta_edited.svg
grep "circle" ${outfolder}/${outtree}_delta_temp.svg >>${outfolder}/${outtree}_delta_edited.svg
tail -1 ${outfolder}/${outtree}_delta_temp.svg >>${outfolder}/${outtree}_delta_edited.svg
perl -pi -e 's/viewBox(.{2})([0-9\.]+) ([0-9\.]+) ([0-9\.]+) ([0-9\.]+)/viewBox.$1.($2-800)." ".($3-800)." ".($4+1600)." ".($5+1600)/ge' ${outfolder}/${outtree}_delta_edited.svg

python convert_svg_to_pdf.py --input ${outfolder}/${outtree}_delta_edited.svg --output ${outfolder}/${outtree}_delta_edited.pdf

