#! /bin/bash

for i in $(cat socs.lst); do
	# get the name of the file
	NAME=$(echo $i | grep -o '^[^.]*')
	SOCCHECK=$(grep -ci "spin_orbit_couple.*true" $NAME.out )
	CMTOHARTREE="0.000000455633"
	if [ $SOCCHECK -ge "1" ]; then
		# How many states per manifold?
		# find the line with the number of roots | keep only the string of numbers at the end of this line
		NSTATES=$(grep "cis_n_roots" $i | sed -e 's/.*[^0-9]\([0-9]\+\)[^0-9]*$/\1/')
		#
		# Singlet to Triplet SOC extraction
		#
		if [ -f $NAME.adiabStoTjmat ]; then # check if adiabatic S to T couplings have been extracted
			echo "Singlet to Triplet couplings already extracted"
		else
			# SOC from singlet GS to triplet excited states
			SEARCHSTRING="Total SOC between the singlet ground state and excited triplet states"
			# find the line that starts the print out, get the next nstates lines, print only the 2nd column * cm-1 to hartree conversion as a row into the matrix as a new file
			awk -v searchstr="$SEARCHSTRING" -v nstates="$NSTATES" -v convert="$CMTOHARTREE" '$0 ~ searchstr {for(l=1; l<=nstates; l++) if (l<nstates) {getline; printf "  %1.5e",$2*convert} else {getline;printf "  %1.5e\n",$2*convert}}' $i > $NAME.temp.adiabStoTjmat
			# Loop through k = S1 to Smax
			for k in $(eval echo "{1..$NSTATES}") ; do
				# Construct the correct search string for state Sk
				SEARCHSTRING="Total SOC between the S$k state and excited triplet states"
				# find the line that starts the print out, get the next nstates lines, print only the 2nd column * cm-1 to hartree conversion as a row into the matrix by appending to end of existing file
				awk -v searchstr="$SEARCHSTRING" -v nstates="$NSTATES" -v convert="$CMTOHARTREE" '$0 ~ searchstr {for(l=1; l<=nstates; l++) if (l<nstates) {getline; printf "  %1.5e",$2*convert} else {getline;printf "  %1.5e\n",$2*convert}}' $i >> $NAME.temp.adiabStoTjmat
			done
			# transpose the whole result
			octave --silent --eval "dlmwrite('$NAME.adiabStoTjmat',dlmread('$NAME.temp.adiabStoTjmat')','  ','precision','% 1.5e')"
			# delete the temp files
			rm -rf $NAME.temp.adiabStoTjmat
		fi
		if [ -f $NAME.adiabToTjmat ]; then # check if adiabatic T to T couplings have been extracted
			echo "Triplet to Triplet couplings already extracted"
		else
			# Touch the temp file so that later concatenation works
			touch $NAME.temp.adiabTtoTjmat
			# Reduce N states by 1 b/c looking for inter-state couplings only
			TSTATES=$(( $NSTATES - 1 ))
			# Loop through the k triplet states
			for k in $(eval echo "{1..$TSTATES}") ; do
				# Construct the correct search string for state Tk
				SEARCHSTRING="Total SOC between the T$k state and excited triplet states"
				# Figure out how many lines to look ahead for this state
				NLINES=$(( $NSTATES - $k ))
				# Figure out width of final matrix in characters with correct padding
				NCHARS=$(( $NSTATES*13 + $NSTATES*2 ))
				# Add rows of zeroes to the correct number of times to the list for a lower triangular matrix
				awk -v searchstr="$SEARCHSTRING" -v nstates="$k" -v convert="$CMTOHARTREE" '$0 ~ searchstr {for(l=1; l<=nstates; l++) {getline; printf "  %1.5e\n",$4*convert}}' $i >> $NAME.temp.adiabTtoTjmat
				# find the search line, print out rows of the next n lines
				awk -v searchstr="$SEARCHSTRING" -v nstates="$NLINES" -v convert="$CMTOHARTREE" '$0 ~ searchstr {for(l=1; l<=nstates; l++) {getline; printf "  %1.5e\n",$2*convert}}' $i >> $NAME.temp.adiabTtoTjmat
			done
			# take the temporary matrix and sort into columns | append a column of zeroes > save to another temp file
			column -c $NCHARS $NAME.temp.adiabTtoTjmat | awk '{print $0, "    0.00000e+00"}' > $NAME.tril.adiabTtoTjmat
			# in octave read in the bottom triangular matrix and add it to the transpose of itself to make the correct symmetric matrix, then save with the correct formatting
			octave --silent --eval "dlmwrite('$NAME.adiabTtoTjmat',dlmread('$NAME.tril.adiabTtoTjmat')+dlmread('$NAME.tril.adiabTtoTjmat')','  ','precision','% 1.5e')"
			# delete the temp files
#			rm -rf $NAME.temp.adiabTtoTjmat
#			rm -rf $NAME.tril.adiabTtoTjmat
		fi
	else
		echo "$NAME.out is not a SOC calculation!"
	fi
done

