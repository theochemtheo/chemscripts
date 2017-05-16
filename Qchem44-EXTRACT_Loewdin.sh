#! /bin/bash

for i in $(cat diabats.lst); do
	# get the name of the file
	NAME=$(echo $i | grep -o '^[^.]*')
        DIABCHECK=$(grep -ci "Entering the .*Localization Code for CIS" $NAME.out )
	if [ $DIABCHECK -ge "1" ]; then
		# How many states in diabatisation?
		NSTATES=$(grep "_cis_numstate" $i | sed -e 's/.*[^0-9]\([0-9]\+\)[^0-9]*$/\1/')
		# Figure out correct character width for the table given the trimming
		NCHARS=$(($NSTATES*14 + $NSTATES*2))
	        # First extract the rotation matrices
	        if [ -f $NAME.rotmat ]; then
			echo "Rotation matrix already extracted from $NAME.out"
	        else
			# sed the correct lines (with optional whitespace) from i | trim the last 14 characters (12 digits + 2 spaces for sign and padding) | columnise the output with a total width of NCHARS | put into output $NAME.rotmat
		        sed -n '/final adiabatic\s*->\s*diabatic/p' $i | grep -o '.\{14\}$' | column -c $NCHARS > $NAME.rotmat
		fi
		if [ -f $NAME.adiabmat ]; then
			echo "Adiabatic energy matrix already extracted from $NAME.out"
	        else
			# grep the correct lines from i | trim the last 14 characters (12 digits + 2 spaces for sign and padding) | columnise the output with a total width of NCHARS | put into output $NAME.rotmat
		        grep "showmatrix adiabatH" $i | grep -o '.\{14\}$' | column -c $NCHARS > $NAME.adiabmat
		fi
		if [ -f $NAME.diabmat ]; then
			echo "Diabatic energy matrix already extracted from $NAME.out"
	        else
			# grep the correct lines from i | trim the last 14 characters (12 digits + 2 spaces for sign and padding) | columnise the output with a total width of NCHARS | put into output $NAME.rotmat
		        grep "showmatrix diabatH" $i | grep -o '.\{14\}$' | column -c $NCHARS > $NAME.diabmat
		fi
	else
		echo "$NAME.out is not a CIS diabatisation calculation!"
	fi
done
