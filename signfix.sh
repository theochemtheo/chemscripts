#! /bin/sh

# Usage signfix.sh [name of mask matrix] [output prefix]
#
# This is a script for fixing arbitrary sign changes in diabatic state matrices due to phase changes in the orbitals, which changes the sign of the CIS amplitudes and thus the phase relations between the excited states
# It utilises the script 'SIGNFIX-DIABMAT.m'

# Store the name of the mask matrix for later
MASKNAME=$1
# Store the output prefix for later
PREFIX=$( echo $2 )

# For the files in sign.lst
for i in $( cat sign.lst ); do
	# find the basename
	# These files must have a name with the format: BLABLA_123_05l.extension
	#                                  1= vib mode number _^1^_^2^ 2= value of dimensionless coordinate Q
	BN=$( echo $i | cut -f 1 -d '.' )
	# Set the output name
	OUTNAME=$( echo "$PREFIX$BN" )
	# Name of the target matrix
	TARMAT=$( echo "$BN.diabmat" )
	# Print progress to terminal
	echo "Reformatting $BN"
	# use the reforming octave script
	octave -qf ~/bin/chemscripts/SIGNFIX-DIABMAT.m $MASKNAME $TARMAT $OUTNAME
done
