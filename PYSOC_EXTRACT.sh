#!/bin/bash

EXSTATES=20
PADDING=8

SOCARRAY=($(awk '{ print $4 }' soc_out.dat))

for SINGLET in $(seq 0 $EXSTATES); do
	SINDEX=$( echo "$SINGLET*$EXSTATES" | bc )
	for TRIPLET in $(seq 1 $EXSTATES); do
		TINDEX=$( echo "$SINDEX+$TRIPLET-1" | bc )
		printf "%${PADDING}.5f\t" "${SOCARRAY[$TINDEX]}"
	done
	printf "\n"
done
