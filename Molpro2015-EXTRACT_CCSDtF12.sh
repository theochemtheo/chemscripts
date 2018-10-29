#!/bin/bash

USAGE="Usage: command calculation.out, where calculation.out is the .out file from a Molpro 2015 CCSD(T)-F12b calculation"


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

# First extract HF values
echo "Hartree Fock energy( or energies)"
grep "  New reference energy" $1
echo "---"
# Then CCSD-F12b correlation energies
echo "CCSD-F12b correlation energy( or energies)"
grep "CCSD-F12b correlation energy" $1
echo "---"
# Triples are printed twice for each calculation, only get the even values
echo "Triples (T) correlation energy( or energies"
grep "Triples (T) contribution" $1 | awk 'NR%2==0'
echo "---"
