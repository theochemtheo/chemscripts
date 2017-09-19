#! /bin/sh

# This script provides an interface between the ABCluster program (PCCP, 2016, 18, 3003)
# and the xtb program (JCTC, 2017, 13, 1989)

# This script is used as follows:
# ./abcluster_gfn-xtb_interface.sh input.xxxxxyz
# where input is in xyz format

# Get name of input file and create variable for the output
BN=$( basename $1 .xxxxxyz )
inname=$( echo "$BN.xyz" )
outname=$( echo "$BN.ooooout" )
logfname=$( echo "$BN.xtb-ouput.out" )
optgeoms=$( echo "$BN.xtb-opt.xyz")
user=$( whoami )
cwd=$( pwd )

# Create a scratch directory for this
SCRATCH=$( echo "/data/$user/gfn-workdir.$$" )

mkdir $SCRATCH

# Copy the input file into the scratch directory
cp $1 $SCRATCH/$inname
cd $SCRATCH

## Run a GFN-xTB optimisation
~/programs/xtb_exe/xtb $SCRATCH/$inname -opt > $SCRATCH/$logfname


# Copy the results into the working directory and clean-up
mv $SCRATCH/xtbopt.coord $cwd/$outname
mv $SCRATCH/$logfname $cwd/
mv $SCRATCH/xtbopt.log $cwd/$optgeoms
rm -rf $SCRATCH
