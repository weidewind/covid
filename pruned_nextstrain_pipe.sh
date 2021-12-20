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
mkdir -p $outfolder
cp $config_file ${outfolder}/$config_file

outtree=tree_from_${from_date}_to_${to_date}
echo "Pruning tree.."
python prune_tree.py --tree $full_tree --dates $leaf_dates_file --from_date $from_date --to_date $to_date --output ${outfolder}/${outtree}.nwk

echo "Drawing pruned tree.."
xvfb-run python draw_nextstrain_tree_sept.py --tree ${outfolder}/${outtree}.nwk --latest_date $to_date --states $leaf_states_file --dates $leaf_dates_file --important_lineages $important_lineages --pangolined $all_pangolined --mode $mode --output ${outfolder}/${outtree}
echo "Editing graphics.."
newr=10
if [ "$mode" = "c" ]; then
	newr=25
fi
echo $newr
sed 's/<g/\n<g/g' ${outfolder}/${outtree}.svg >${outfolder}/${outtree}_temp.svg
sed -i "s/r=\"2\"/r=\"$newr\"/g" ${outfolder}/${outtree}_temp.svg
grep -v "circle" ${outfolder}/${outtree}_temp.svg | head -n -1 >${outfolder}/${outtree}_edited.svg
grep "circle" ${outfolder}/${outtree}_temp.svg >>${outfolder}/${outtree}_edited.svg
tail -1 ${outfolder}/${outtree}_temp.svg >>${outfolder}/${outtree}_edited.svg
perl -pi -e 's/viewBox(.{2})([0-9\.]+) ([0-9\.]+) ([0-9\.]+) ([0-9\.]+)/viewBox.$1.($2-800)." ".($3-800)." ".($4+1600)." ".($5+1600)/ge' ${outfolder}/${outtree}_edited.svg

python convert_svg_to_pdf.py --input ${outfolder}/${outtree}_edited.svg --output ${outfolder}/${outtree}_edited.pdf

