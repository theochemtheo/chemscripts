#! /bin/bash

USAGE="Usage: command input, where input is the name of a (non-spin-flip) Qchem 5.0 calculation including the CIS/TDA derivative couplings"

# For details regarding these calculations, please see J. Chem. Phys, 2011, 135, 234105
# This script will create three matrices for each pair of states
#
# Theo Keane, 19 March 2017

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

# check this is a Qchem 5.0 calculation
if grep -q "Q-Chem 5.0, Q-Chem, Inc., Pleasanton, CA (2017)" $1; then
# check this includes derivative couplings
if grep -i -q "cis_der_numstate" $1; then
# Begin extracting the states
	echo "Extracting the derivative coupling vectors from $1"
	# Start defining important variables
	NAME=$( echo $1 | grep -o '^[^.]*' )
	# Figure out number of atoms
	NATOMS=$( grep -B2 "Nuclear Repulsion Energy" $1 | head -n1 | awk '{print $1;}' )
	# Since only one spin manifold can be calculated at a time, check multiplicity of excited states via triplets
	TRIPSCHECK=$( grep -i "cis_triplets" $1 | awk 'NF{ print $NF }' )
	if [[ "$TRIPSCHECK" = *"true"* ]]; then
		EXMANIFOLD=T
	else
		EXMANIFOLD=S
	fi
	# find which derivative couplings are involved
	DCSTATES=$( awk '/derivative_coupling/{getline; getline; print}' $1 )
	# loop over the pairs of states
	set -- $DCSTATES
	for STATEI; do
		# Define the multiplicity of state i
		STATEIMULT=$EXMANIFOLD
		# Check if GS is involved
		if [ $STATEI = 0 ]; then
			# Find multiplicity of GS
			GSMULTCHECK=$( awk '/\$molecule/{ getline; print $NF }' $NAME.out )
			# Amend multiplicty variable for the multiplicity of the GS
			if [[ "$GSMULTCHECK" = *"1"* ]]; then
				STATEIMULT=S
			else
				STATEIMULT=T
			fi
		fi
		shift
		for STATEJ; do
			# Stem of the name for the output
			OUTPUTSTEM=$( echo "$NAME.$STATEIMULT$STATEI-to-$EXMANIFOLD$STATEJ" | tr -cd '[:print:]' )
			# the string to be searched for for this pair of states
			SEARCHSTR=$( printf "between states %s and %s" "$STATEI" "$STATEJ" )
			echo $SEARCHSTR
			# prepare to search for DC without ETF
			# line to start reading out
			READOUTSTART=10
			# line to stop reading out
			READOUTEND=$(( $READOUTSTART + $NATOMS ))
			awk -v s="$SEARCHSTR" -v n="$READOUTSTART" -v m="$READOUTEND" '$0 ~ s {for(l=1; l<=m; l++) if (l<=n) {getline;} else {getline; printf "% 1.6f\t% 1.6f\t% 1.6f\n", $2, $3, $4} }' $NAME.out > $OUTPUTSTEM.DCnoETF.mat
			# Adjust params to for GD
			READOUTSTART=$(( $READOUTEND + 5 ))
			READOUTEND=$(( READOUTSTART + $NATOMS ))
			awk -v s="$SEARCHSTR" -v n="$READOUTSTART" -v m="$READOUTEND" '$0 ~ s {for(l=1; l<=m; l++) if (l<=n) {getline;} else {getline; printf "% 1.6f\t% 1.6f\t% 1.6f\n", $2, $3, $4} }' $NAME.out > $OUTPUTSTEM.GD.mat
			# Adjust params for DC + ETF
			READOUTSTART=$(( $READOUTEND + 5 ))
			READOUTEND=$(( READOUTSTART + $NATOMS ))
			awk -v s="$SEARCHSTR" -v n="$READOUTSTART" -v m="$READOUTEND" '$0 ~ s {for(l=1; l<=m; l++) if (l<=n) {getline;} else {getline; printf "% 1.6f\t% 1.6f\t% 1.6f\n", $2, $3, $4} }' $NAME.out >  $OUTPUTSTEM.DCwithETF.mat
		done
	done
else
# does not containt cis_der_numstate
	echo "$1 does not include derivative couplings"
fi
else
# is not a Qchem 5.0 calculation
	echo "$1 is not a Qchem 5.0 calculation"
fi

