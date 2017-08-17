#!/usr/bin/python3

import argparse
import re
import numpy as np
from numpy import linalg as LA


# Ensure line width is wide enough for matrices
np.set_printoptions(linewidth=160)


# Function for obtaining the absolute energy of state i
def getEXE(i):
    EXEi = np.float64(re.findall('(?<=Total energy for state{0:4}:).*'.format(i), open(args["file"]).read())[0])
    return EXEi


# Set up argparse
parser = argparse.ArgumentParser(description="This script is for converting Qchem 4.4 derivative couplings (d_JK) to non-adiabatic couplings (h_JK)")
parser.add_argument("-i", dest="file", metavar="input", help="the Qchem 4.4 output file", required=True)
parser.add_argument("-J", dest="Sone", metavar="<state J>", help="the first state, J, of h_JK", type=int, required=True)
parser.add_argument("-K", dest="Stwo", metavar="<state K>", help="the second state, K, of h_JK", type=int, required=True)
parser.add_argument("-v", dest="show", metavar="<True/False>", help="should the calculated vector h_JK be shown? (must be capitalised)", type=bool, default=False)
parser.add_argument("-o", dest="outp", metavar="output", help="output destination of the vector if desired", required=False)
args = vars(parser.parse_args())

# Global variables
# # hartrees to eV
# autoeV = 27.2114

# Figure out the number of excited states from keywords
ndim = int(re.findall('(?<=cis_n_roots).*[0-9]+', open(args["file"]).read())[0])
# Are states J and K even calculated here?
if args["Sone"] > ndim:
    print("State {} is not one of the roots in this calculation".format(args["Sone"]))
    exit()
elif args["Stwo"] > ndim:
    print("State {} is not one of the roots in this calculation".format(args["Stwo"]))
    exit()

# Figure out the number of atoms
dumbo = str(re.findall('((.*\n){2}) Nuclear Repulsion Energy', open(args["file"]).read())[0])
natoms = int(re.findall('[0-9]+', dumbo)[0])

# J must be a lower state than K for the search string to work
if args["Sone"] < args["Stwo"]:
    J = args["Sone"]
    K = args["Stwo"]
else:
    J = args["Stwo"]
    K = args["Sone"]

# # Extract the ground state energy
# GSE = np.float64(re.findall('(?<=Total energy in the final basis set =).*', open(args["file"]).read())[0])

# EXenergies = []
# for i in range(1, ndim + 1):
#     EXenergies.append(getEXE(i))

# print("The file is {}".format(args["file"]))
# print("The dimensionality is {}".format(ndim))
# print("{}".format(EXenergies))

# Create the search string for this pair of states
searchstr = "between states {} and {}".format(J, K)
# Calculate difference in energy between J and K
EdiffJK = getEXE(K) - getEXE(J)
# prepare empty normal list that will contain the extracted part of the file
allcouplingtext = []
# loop through the lines of the file until the search string is found
with open(args["file"]) as fh:
    for line in fh:
        if searchstr in line:
            # when the search string is found, copy all lines that are associated with the output into the list
            for i in range(natoms * 3 + 21):
                # remove the newline characters and keep only the first 4 columns (NAtom, X, Y, Z), as a list of strings
                allcouplingtext.append((next(fh).rstrip()).split()[0:4])
# Once the relevant part of the file has been extracted, clean it up into a 4 column (NAtom, X, Y, Z) np array of floats and put into individual arrays
# DC no ETF: the derivative coupling vector, d_JK = < Psi_J | d/dR | Psi_K > = h_JK / (E_J - E_K). This is unused at this time
# DCnoETF = np.array(allcouplingtext[10:10 + natoms]).astype('f')
# NACV: the non-adiabatic coupling vector, h_JK, = < Psi_J | dH/dR | Psi_K >. This is usually quite low precision because the output formatting is '-0.6f'
NACV = np.array(allcouplingtext[15 + natoms: 15 + natoms * 2]).astype('f')
# DCwithETF: the derivative coupling vector, d_JK, including the electron-ranslation factor correction. See J.Chem.Phys., 2011, 135, 234105.
DCwithETF = np.array(allcouplingtext[20 + natoms * 2: 20 + natoms * 3]).astype('f')

# Re-calculate the NACV using the DCwithETF and delta E to recover a higher precision number
calcNACV = np.multiply(DCwithETF[:, 1:4], EdiffJK)

# Output the magnitudes to the screen
print("The NACV magnitude obtained directly from the file is:             {:+.4e}".format(LA.norm(NACV[:, 1:4], 'fro')))
print("The NACV magnitude as derived from the DC with ETF and delta E is: {:+.4e}".format(LA.norm(calcNACV, 'fro')))

# If the optional condition is set, print the whole NACV to the screen
if args["show"] is True:
    print("The NACV derived from the DC with ETF and delta E is:")
    print("{}".format(calcNACV))

# If the optional condition is set, print the whole NACV to the desired file
if args["outp"] is not None:
    print("The calculated NACV will be printed to file: {}".format(args["outp"]))
    np.savetxt(args["outp"], calcNACV, fmt="%+.6e", delimiter='\t')
