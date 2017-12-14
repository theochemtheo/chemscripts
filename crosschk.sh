#! /bin/bash

# Usage crosschk.sh [output prefix] [SCF energy in a.u. which will be relative 0]
# This is a script for checking if there has been a diabatic state crossing over a set of input diabatic state matrices
# It utilises the Octave script 'CROSSOVER-CHECK.m'

# Store the output prefix for later
PREFIX=$1

# For the files in perm.lst
for i in $( cat perm.lst ); do
	# find the basename
	# These files must have a name with the format: BLABLA_123_05l.diabmat
	#                                  1= vib mode number _^1^_^2^ 2= value of dimensionless coordinate Q
	BN=$( basename $i .diabmat )
	# find the l or r suffix
	SUFFIX=$( echo $BN | awk '$0=$NF' FS= )
	# Name of the output
	OUTMAT=$( echo "$1_$BN.diabmat" )
	# Construct the name of step |Q|-2
	ONAME=$( echo $BN | awk -v p="$PREFIX" -v s="$SUFFIX" -F_ 'BEGIN{OFS="_";} {j=$(NF)-2;$NF=""; if(j <= 0) s="q"; if(j <= 0) j="0"; printf "%s_%s%02d%s.diabmat",p,$0,j,s}' )
	# Construct the name of step |Q|-1
	NNAME=$( echo $BN | awk -v p="$PREFIX" -v s="$SUFFIX" -F_ 'BEGIN{OFS="_";} {j=$(NF)-1;$NF=""; if(j <= 0) s="q"; if(j <= 0) j="0"; printf "%s_%s%02d%s.diabmat",p,$0,j,s}' )
	# If MNAME = NNAME = ONAME, i.e. if all are BLABLA_123_00q.diabmat, no manipulation should be done, instead just copy to the desired output matrix and skip the rest of this loop
	if [ "$ONAME" == "$NNAME" ] && [ "$NNAME" == "$OUTMAT" ]; then
		echo "$i cannot be checked for a crossover, copying $i to $OUTMAT"
		cp $i $OUTMAT
		continue
	fi
	# Relative energy for this point
	EREL=$( awk -v s="$2" '/ion met/ {printf "% 1.10f",$2-s}' $BN.out )
	# Print progress to terminal
	echo "Reformatting $BN"
	# use the reforming octave script
	octave -qf ~/bin/chemscripts/CROSSOVER-CHECK+EREL.m $ONAME $NNAME $i $OUTMAT $EREL
done
