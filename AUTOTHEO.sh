#!/bin/bash

if [[ -z "${1}" ]]; then
  echo "You must specify which atom is the metal"
  exit 1
fi

if [[ ! -v $CHEMSCRIPTS ]]; then
  export CHEMSCRIPTS=/home/ch1tk/bin/chemscripts
fi

export EXPECTSCRIPTS=${CHEMSCRIPTS}/expect_scripts

# First prepare sing_a and trip_a
if [[ -f ciss_a ]] && [[ ! -f sing_a ]]; then
  echo "Creating softlink to ciss_a"
  ln -s ciss_a sing_a
fi

if [[ -f cist_a ]] && [[ ! -f trip_a ]]; then
  echo "Creating softlink to cist_a"
  ln -s cist_a trip_a
fi

if [[ ! -f molden.input ]]; then
  echo "No molden.input found, running tm2molden"
  ${EXPECTSCRIPTS}/tm2molden.exp
fi

echo "Running theoinp"
${EXPECTSCRIPTS}/theoinp_escf.exp ${1}

echo "Running analyze_tden.py"
nohup analyze_tden.py >/dev/null 2>&1
