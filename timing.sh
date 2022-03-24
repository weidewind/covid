#!/usr/bin/bash


function timing {
		input=$1
		# Internal nodes dating
		# trying to get dates of foreign and russian strains, which are closest to entry node (entry node dating for the poor)
		echo "Foreign strains timing.."
		start=`date +%s`
		#split file with entries into chuncks
		entry_num=$(grep -c -v "^$" $input)
		echo "there are $entry_num lines in $input file"

		chunk_size=$((entry_num/threads))
		filenum=$((threads+1))
		echo "chunk size is $chunk_size lines"
		for i in $( seq 1 $threads )
			do
			echo "head -$((i*chunk_size)) $input | tail -$chunk_size >$input.part${i}"
			head -$((i*chunk_size)) $input | tail -$chunk_size >$input.part${i}
			done
		tail -$((entry_num-threads*chunk_size)) $input >$input.part${filenum}

		for i in $( seq 1 $((threads+1)) )
		do
			( python -u timing.py --tree ${all_states_file}.newick.broken --dates $leaf_dates_file --countries $leaf_states_file  --entry_nodes_file $input.part${i} --output ${input}.part${i} >${input}.part${i}.dates_stats.log ) &
		done
		wait

		cat ${input}.part1.dates_stats >${input}.dates_stats
		for i in $( seq 2 $((threads+1)) )
			do
			tail -n +2 ${input}.part${i}.dates_stats >>${input}.dates_stats
			#unlink ${input}.part${i}.dates_stats
			#unlink ${input}.part${i}.dates_stats.log
			done

		end=`date +%s`
		runtime=$( echo "$end - $start" | bc -l )
		echo "Entry timing performed in "$((runtime/60))" min "
}
