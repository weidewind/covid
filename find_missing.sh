#!/usr/bin/bash
while getopts i:o:t: option
do
		case "${option}" in
			i) infolder=${OPTARG};; # folder with fasta
			o) outfolder=${OPTARG};; # folder with .pangout files
			t) tree=${OPTARG};;
		esac
done

find $infolder -regex ".*\.fa[a-z]*" -print0 | while read -d $'\0' file
do
   fname=$(basename $file)
   if [ ! -f "${outfolder}/${fname}.pangout" ]; then
      strain=$(echo $fname |cut -d'.' -f1)
      if grep -q "[,(]${strain}" $tree; then
         echo $fname
      fi
   fi
done