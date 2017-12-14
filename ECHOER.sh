#! /bin/bash

# Usage: echoer.sh [total number of lines per iteration] [output name]

# delete previous output
rm -f $2

# Initialize counter
counter=1

# Loop through the files that match the pattern
for i in $(cat elist.lst); do
	# For the first loop, check how long the matching files are
	while [ $counter == 1 ]; do
		# Number of lines for the matching file
		NLINES=$( wc -l < $i )
		# Calculate number of newlines required
		REPLINES=$(($1-$NLINES))
		# Expand this
		REPS=$(eval echo {1..$REPLINES})
		# Increment the counter
		((counter++))
	done
	cat $i >> $2
	printf '\n%.0s' $REPS >> $2
	echo $i
done
