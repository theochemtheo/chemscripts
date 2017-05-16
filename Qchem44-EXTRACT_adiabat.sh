#!/bin/sh

USAGE="Usage: This script requires a file 'diabats.lst' which must contain the names of Qchem 4.4 calculations including CIS/TDA excitations."

# Theo Keane, 31 March 2017

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

for i in $(cat diabats.lst); do
	# get the name of the file
	NAME=$( echo $i | grep -o '^[^.]*')
	# check this is a Qchem 4.4 calculation
	QCHECK=$(grep -ci "Q-Chem 4.4, Q-Chem, Inc., Pleasanton, CA (2016)" $NAME.out )
	if [ $QCHECK -ge "1" ]; then
	# check this includes CIS excitations
	CISCHECK=$(\grep -c "Excitation Energies" $NAME.out )
	if [ $CISCHECK -ge "1" ]; then
	# Begin extracting the states
		# Name of energies output
		ENOUTPUT=$( echo "$NAME.adiabex" )
		# Name of dipole output
		MUOUTPUT=$( echo "$NAME.dipoles" )
		# check if excitations have already been extracted
		if [ -f $ENOUTPUT ]; then
			echo "Adiabatic excitations already extracted from $i!"
		else
		echo "Extracting the adiabatic excitation energies from $i"
		# find how many adiabatic states there are
		NSTATES=$( awk '/cis_n_roots/{print $NF}' $i )
		# GS energy in Hartrees
		GSE=$( awk '/Convergence criterion met/{print $2}' $i )
		# Touch the output to ensure that future concatenations work OK
		touch $ENOUTPUT
		for ((j=1;j<=$NSTATES;j++)); do
			# Set up the search string for this state
			printf -v SEARCHSTR 'Total energy for state%4s' $j
			awk -v s="$SEARCHSTR" -v GSE="$GSE" '$0 ~ s {printf "% 1.12f\n", $6-GSE}' $NAME.out >> $ENOUTPUT
		done
		fi
		if [ -f $MUOUTPUT ]; then
			echo "Dipoles already extracted from $i!"
		else
		echo "Extracting dipoles from $i"
		# find how many adiabatic states there are
		NSTATES=$( awk '/cis_n_roots/{print $NF}' $i )
		# Touch the output to ensure that future concatenations work OK
		touch $MUOUTPUT
		for ((j=1;j<=$NSTATES;j++)); do
			# Set up the search string for this state
			printf -v SEARCHSTR 'Excited-State Multipoles, State%3s' $j
			awk -v s="$SEARCHSTR" '$0 ~ s {getline;getline;getline;getline;getline;getline;printf "% 1.4f\n", $2}' $NAME.out >> $MUOUTPUT
		done
		fi
	else
	# does not contain cis_n_roots
		echo "$i does not include CIS(RPA)/TDA-(TD-)DFT excitations"
	fi
	else
	# is not a Qchem 4.4 calculation
		echo "$i is not a Qchem 4.4 calculation"
	fi
done
