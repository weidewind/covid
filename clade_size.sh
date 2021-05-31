#!/usr/bin/bash

mapfile -t lineages < <( cut -f2 pangoline.all | grep -E "B\.1\.1\.[0-9]+$" | sort | uniq )
for lin in ${lineages[*]}
	do
	count=$( cut -f2 pangoline.all | grep -c -E "$lin$|$lin\." )
	echo "$lin $count"
	done
