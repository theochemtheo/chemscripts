#! /bin/bash

USAGE="Usage: command input, where input is a Gaussian 09 frequency calculation output file, for a calculation including the keyword freq=hpmodes"

# For more information about these types of modes, search for "Vibrational Analysis in Gaussian". You should be able to find it on the Gaussian website.
# Alternatively, read "Molecular Vibrations: The Theory of Infrared and Raman Vibrational Spectra", by Wilson, Decius and Cross. (Dover, 1955).
#
# The calculation must be run with the freq=hpmodes keyword for this to work.
#
# The matrices are deposited into files named [log-file-name].m[n].dat, where [n] is the number of the mode. If you have a big molecule, this will mean a lot of files are produced.
#
# Theo Keane, 15 March 2017

# Make sure a file has been specified
if [ $# -eq 0 ]; then
	echo "$USAGE"
	exit 1;
fi

# Check for a -h flag
while getopts ":h" optname
do
	case "$optname" in
		"h")
		echo $USAGE
		exit 0;
		;;
		"?")
		echo "There are no options for this script"
		echo $USAGE
		exit 0;
		;;
		*)
		echo "Unknown error while processing options"
		echo $USAGE
		exit 0;
		;;
	esac
done
# End of checking for flags

# Check that this is a Gaussian 09 calculation
if grep -q "This is part of the Gaussian(R) 09 program." $1; then
# If it is a Gaussian 09 calculation, check that the file includes the keyword freq=hpmodes
if grep -q 'freq.*hpmodes' $1; then
# If it does contain freq=hpmodes, continue
echo "Extracting the high-precision, (Frobenius) normalized, mass-weighted cartesian displacement matrices for the normal modes"
# Extract the name, number of degrees of freedom, number of atoms and the number of lines to be searched
NAME=$(echo $1 | grep -o '^[^.]*')
NATOMS=$( grep -o -m 1 "NAtoms.*NActive" $1 | sed -e 's/.*[^0-9]\([0-9]\+\)[^0-9]*$/\1/' )
NLINES=$(( 3 * $NATOMS ))
DOF=$(( 3 * $NATOMS - 6 ))
# Loop over the DOFs
for i in $(seq 1 $DOF);do
	# Which index of search term should be used for this mode?
	SINDEX=$(( $i / 5 + 1 ))
	# Modulo 5 the DOF index
	MOD=$((  $i % 5 ))
	# If index is a multiple of 5 then correct MOD and SINDEX, plus announce progress
	if [ $MOD = 0 ]; then
		# Set the MOD variable to 0 + 5 = 5
		MOD=$(( $MOD + 5 ))
		SINDEX=$(( ($i - 1) / 5  + 1))
		echo "Extracting mode $i"
	fi
	# Set which column is to be searched
	COLUMN=$(( 3 + $MOD ))
	# Search for SINDEX'th entry of matching preceeding line, then for the next NLINES take the correct COLUMN and print with tabs unless it's a z entry in which case print with tab and newline
	awk -v m="$NLINES" -v s="$SINDEX" -v c="$COLUMN" '$0 ~ "Coord Atom Element:" && !--s {for(l=1; l<=m; l++) if (l<=m) if (l%3<1) {getline; printf "\t% 1.5f\n", $c} else {getline; printf "\t% 1.5f", $c} }' $1 > $NAME.m$i.dat
done
echo "Extraction Complete! If your molecule is linear, the final mode was skipped - make sure to extract it yourself manually!"
else
# File does not include freq=hpmodes
echo "The calculation must include freq=hpmodes"
fi
else
# File is not a Gaussian 09 calculation
echo "The file $1 is not a Gaussian 09 calculation"
fi
