#!/usr/bin/python3
"""
This script is intended to provide examples of ways of getting excitations out of the output from a Turbomole 7.2 escf calculation

These are intended to emulate the types of array you get from cclib
"""

import argparse
import numpy as np
import re


parser = argparse.ArgumentParser()
parser.add_argument("escffile", help="name of escf output file")
args = parser.parse_args()

excitations = np.loadtxt(re.findall('(?<=Excitation energy \/ cm\^\(-1\)\:).*', open(args.escffile).read()))
velrep = np.loadtxt(re.findall('(?<=Oscillator strength\:\n\n    velocity representation\:).*', open(args.escffile).read()))
lenrep = np.loadtxt(re.findall('(?<=length representation\:).*', open('escf.out').read()))[0::2]
mixedrep = np.loadtxt(re.findall('(?<=mixed representation\:).*', open(args.escffile).read()))

print("{}".format(excitations))
print("{}".format(velrep))
print("{}".format(lenrep))
print("{}".format(mixedrep))
