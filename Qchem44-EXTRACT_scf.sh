#! /bin/bash

for i in $(cat diabats.lst); do
	# get the name of the file
	NAME=$(echo $i | grep -o '^[^.]*')
        DIABCHECK=$(grep -ci "Entering the .*Localization Code for CIS" $NAME.out )
	if [ $DIABCHECK -ge "1" ]; then
		grep "criterion met" $i | awk '{print $2}' > $NAME.scf
	else
		echo "$NAME.out is not a CIS diabatisation calculation!"
	fi
done
