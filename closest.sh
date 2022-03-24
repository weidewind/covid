#!/usr/bin/bash


function closest {
		input=$1
		# echo "Searching for closest ingroup and outgroup strains for nodes listed in $input .."
		start=`date +%s`
		#split file with entries into chuncks
		entry_num=$(grep -c -v "^$" $input)
		# echo "there are $entry_num lines in $input file"

		chunk_size=$((entry_num/threads))
		filenum=$((threads+1))
		# echo "chunk size is $chunk_size lines"
		for i in $( seq 1 $threads )
			do
			# echo "head -$((i*chunk_size)) $input | tail -$chunk_size >$input.part${i}"
			head -$((i*chunk_size)) $input | tail -$chunk_size >$input.part${i}
			done
		tail -$((entry_num-threads*chunk_size)) $input >$input.part${filenum}

		for i in $( seq 1 $((threads+1)) )
		do
			( python -u closest.py --tree ${all_states_file}.newick.broken --nodes_file $input.part${i} --output ${input}.inout.part${i} >${input}.inout.part${i}.dates_stats.log ) &
		done
		wait

		cat ${input}.inout.part1.dates_stats >${input}.inout.dates_stats
		for i in $( seq 2 $((threads+1)) )
			do
			tail -n +2 ${input}.inout.part${i}.dates_stats >>${input}.inout.dates_stats
			#unlink ${input}.part${i}.dates_stats
			#unlink ${input}.part${i}.dates_stats.log
			done

		end=`date +%s`
		runtime=$( echo "$end - $start" | bc -l )
		# echo "Search performed in "$((runtime/60))" min "
		# echo "\n"
		echo ${input}.inout.dates_stats
}

function closest_for_random {
		output=$1
		entry_num=$2
		# echo "Searching for closest ingroup and outgroup strains for nodes listed in $input .."
		start=`date +%s`
		#split file with entries into chuncks
		# echo "there are $entry_num lines in $input file"
		chunk_size=$((entry_num/threads))

		for i in $( seq 1 $((threads+1)) )
		do
			( python -u closest.py --tree ${all_states_file}.newick.broken --shift --random $chunk_size --output ${output}random.inout.part${i} >${output}random.inout.part${i}.dates_stats.log ) &
		done
		wait

		cat ${output}random.inout.part1.dates_stats >${output}random.inout.dates_stats
		for i in $( seq 2 $((threads+1)) )
			do
			tail -n +2 ${output}random.inout.part${i}.dates_stats >>${output}random.inout.dates_stats
			#unlink random.part${i}.dates_stats
			#unlink random.part${i}.dates_stats.log
			done

		end=`date +%s`
		runtime=$( echo "$end - $start" | bc -l )
		# echo "Search performed in "$((runtime/60))" min "
		# echo "\n"
		echo ${output}random.inout.dates_stats
}