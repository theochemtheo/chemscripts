#!/bin/bash

# check_dens_ana.sh [name of dens_ana.in] [number of atoms]

DENSANA=$1
NATM=$2

# First check number of atoms
COMMACOUNT=$( grep "at_lists" ${DENSANA} | grep -o "," | wc -l )
ATOMCOUNT=$( echo "$COMMACOUNT + 1" | bc )

if [ $ATOMCOUNT -eq $NATM ]; then
  echo "Atom count correct"
else
  echo "Atom count in $DENSANA is $ATOMCOUNT, not $NATM"
fi

echo "Checking for duplicates and missing atoms"
errorcount=0
for i in $(seq 1 ${NATM}); do
  NATMCOUNT=$( grep -Eo "(\[${i}\],|\[${i},|, ${i}\]|, ${i},|,${i},|,${i}\])" $DENSANA | wc -l )
  if [ $NATMCOUNT -gt 1 ]; then
    errorcount+=1
    echo "Warning, atom label ${i} appears at least ${NATMCOUNT} times!"
  elif [ $NATMCOUNT -eq 0 ]; then
    errorcount+=1
    echo "Warning, atom label ${i} is missing!"
  fi
done
if [ $errorcount -eq 0 ]; then
  echo "No duplicates or missing atoms found"
fi
