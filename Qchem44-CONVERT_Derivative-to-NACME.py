#!/usr/bin/python3

import argparse
import re
import numpy as np


# Ensure line width is wide enough for matrices
np.set_printoptions(linewidth=160)


def getEXE(i):
    EXEi = np.float64(re.findall('(?<=Total energy for state{0:4}:).*'.format(i), open(args["file"]).read())[0])
    return EXEi


parser = argparse.ArgumentParser(description="This script is for converting Qchem 4.4 derivative couplings (d_JK) to non-adiabatic couplings (h_JK)")
parser.add_argument("-f", dest="file", metavar="file", help="the Qchem 4.4 output file", required=True)
args = vars(parser.parse_args())

# Global variables
# # hartrees to eV
# autoeV = 27.2114

# Figure out the number of excited states from keywords
ndim = int(re.findall('(?<=cis_n_roots).*[0-9]+', open(args["file"]).read())[0])
# Figure out the number of atoms
dumbo = str(re.findall('((.*\n){2}) Nuclear Repulsion Energy', open(args["file"]).read())[0])
natoms = int(re.findall('[0-9]+', dumbo)[0])


# # Extract the ground state energy
# GSE = np.float64(re.findall('(?<=Total energy in the final basis set =).*', open(args["file"]).read())[0])

EXenergies = []
for i in range(1, ndim + 1):
    EXenergies.append(getEXE(i))

# print("The file is {}".format(args["file"]))
# print("The dimensionality is {}".format(ndim))
print("{}".format(EXenergies))

searchstr = "between states 1 and 2"
allcouplingtext = np.empty([0, natoms * 3 + 14], dtype='U67')
with open(args["file"]) as fh:
    for line in fh:
        if searchstr in line:
            for i in range(natoms * 3 + 14):
                # allcouplingtext.append(next(fh))
                allcouplingtext = np.insert(allcouplingtext, i, next(fh))

print("{}".format(allcouplingtext))
