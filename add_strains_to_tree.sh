#!/usr/bin/bash

# Take global gisaid phylogeny and gisaid msa, 
# ?? exclude any sequences not present in the tree ??
# add russian strains to gisaid msa,
# remove from the tree any strains that are missing in the resulting alignment (if there are any)
# add to the tree sequences that are missing from the tree but are present in the alignment
# mostly copied from https://github.com/roblanf/sarscov2phylo/blob/master/scripts/global_tree_gisaid_start_tree.sh
helpFunction()
{
  # echo "Make an ML phylogeny from a large FASTA file of GISAID sequences"
  # echo "Usage: $0 -i GISAID_fasta -o output_filename -s start_tree -t threads "
  # echo "    -i Full path to unaligned fasta file of SARS-CoV-2 sequences from GISAID"
  # echo "    -p Filpath to previous iteration to be updated (must contain an ft_SH.tree file and an excluded_sequences.tsv file)"
  # echo "    -t number of threads to use"
   exit 1 # Exit script after printing help
}

while getopts "t:f:g:n:" opt
do
   case "$opt" in
      t ) inputtree="$OPTARG" ;;
	  f ) inputfasta="$OPTARG" ;;
      g ) globmsa="$OPTARG" ;;
      n ) threads="$OPTARG" ;;
      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done


# Print helpFunction in case parameters are empty
if [ -z "$inputfasta" ] || [ -z "$inputtree" ] || [ -z "$globmsa" ] || [ -z "$threads" ]
then
   echo "Some or all of the parameters are empty";
   helpFunction
fi

DIR="$(cd "$(dirname "$0")" && pwd)"

inputdir=$(dirname $inputfasta)
cd $inputdir
mkdir tmp


# 1. Take the alignment and the tree, remove the irrelevant sequences with R

echo ""
echo "Removing unused sequences from input tree"
echo ""
grep ">" $globmsa | cut -c 2- > tmp/alignment_names.txt
Rscript $DIR/clean_tree.R $inputtree tmp/alignment_names.txt


# 2. Todo add our russian strains to the global alignment

echo "Aligning all new seuqences to the global alignment profile"

echo "Cutting reference sequence into tmp/reference.fa"
# cut WH1 from the global alignment and use it as a reference (yes, align to reference and not to alignment profile)
# https://mafft.cbrc.jp/alignment/software/closelyrelatedviralgenomes.html

#EPI_ISL_406798,hCoV-19/Wuhan/WH01/2019,2019-12-26
refid="EPI_ISL_406798"
echo $refid >reference.txt
faOneRecord $globmsa $refid >reference.fa

# so we can get at it as a global variable
export REFERENCE_ALN=reference.fa
export REFERENCE_TXT=reference.txt

echo ""
echo "Splitting new sequences into ${threads} chunks"
faSplit sequence $inputfasta $threads msa_chunk

echo ""
echo "Profile aligning each chunk to the target alignment"
echo "This can take a while, be patient"
echo ""

profile_align()
{

	seqfile=$1
	alfile=$seqfile"_profile_aligned.fa"
	#stofile=$alfile".sto" # dunno why this is here, never used in the original script
    final=$seqfile"_ind_aligned.fa"

	# probably use addfragments instead? should be faster. But it's not, in fact, it is muuuch more slow
	# https://mafft.cbrc.jp/alignment/software/closelyrelatedviralgenomes.html
	# mafft --thread 1 --addfragments fragments --thread -1 existing_alignment > output

	
	START=$(date +%s.%N)
	mafft --thread 1 --quiet --keeplength --add $seqfile "$REFERENCE_ALN" > $alfile
	END=$(date +%s.%N)
	DIFF=$(echo "$END - $START" | bc)
	echo "mafft --add to a reference seq done in ${DIFF} sec" >>mafft_add.timelog
	
	#name=$(grep ">" $seqfile | tr -d ">")
	#echo "$name" | faSomeRecords $alfile /dev/stdin $final
	
	
	faSomeRecords -exclude $alfile $REFERENCE_TXT $final

	rm $seqfile
	rm $alfile

}

export -f profile_align
touch mafft.timelog

inputdir=$(dirname $inputfasta)
ls $inputdir | grep msa_chunk | parallel -j $threads --bar "profile_align {}" > /dev/null

cp $globmsa ${globmsa}plus
# join it all together and clean up
# note that here we *APPEND* to the global alignment, which allows us to add small batches of new stuff whenever we like
ls *_ind_aligned.fa | xargs cat >> ${globmsa}plus
find $inputdir -maxdepth 1 -name "*_ind_aligned.fa" -delete


# 3. Add the new seuqences with IQ-TREE
echo ""
echo "Adding new sequences to input tree with IQ-TREE"
echo ""

# this just adds the new sequences with parsimony
# benchmarking shows that 1 thread is optimal
iqtree2 -seed 1729 -s ${globmsa}plus -g input_tree_cleaned.tree -n 0 -m JC -fixbr -nt 1 --suppress-zero-distance --suppress-list-of-sequences --suppress-duplicate-sequence -pre iqtree_seqsadded_mp
echo ""
echo "Optimising tree with fasttree MP"
echo ""
# we have to do some contortions to set the optimal number of threads for fasttree, which is 3 (see fasttreeOMP.md)
env > old_env.txt
old_threads=$(grep -hoP '^OMP_NUM_THREADS=\K\d+' old_env.txt)
rm old_env.txt
export OMP_NUM_THREADS=3
FastTreeMP -nt -gamma -nni 0 -spr 2 -sprlength 1000 -boot 100 -log fasttree.log -intree iqtree_seqsadded_mp.treefile ${globmsa}plus > $globmsa'_ft_SH.tree'
if [ -n "$old_threads" ]; then
    export OMP_NUM_THREADS=$old_threads
else
    unset OMP_NUM_THREADS
fi

#echo ""
#echo "Cleaning tree with treeshrink"
#echo ""
#run_treeshrink.py -t $globmsa'_ft_SH.tree' -q 0.05 -c -o treeshrink_SH

echo ""
echo "Re-rooting tree on hCoV-19/Wuhan/WH04/2020|EPI_ISL_406801|2020-01-05"
echo "see https://www.biorxiv.org/content/10.1101/2020.04.17.046086v1"
echo ""
#nw_reroot 'treeshrink_SH/'$globmsa'_ft_SH_0.05.tree' "'EPI_ISL_406801'" > ft_SH.tree
python $DIR/reroot.py --tree $globmsa'_ft_SH.tree' --root $refid --output "expanded_tree.newick"


#sed -i.bak "s/'//g" ft_SH.tree
#rm ft_SH.tree.bak