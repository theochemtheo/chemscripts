#!/usr/bin/python3

'''
This script uses parseQCHEM to extract non-adiabatic coupling vectors from an output file and prints data to the screen
'''


import argparse
import numpy as np
import parseQCHEM as pQ


parser = argparse.ArgumentParser(description="This script prints non-adiabatic couplings from a Q-Chem 5.0 file")
parser.add_argument("-i", dest="file", metavar="file", help="input Q-Chem 5.0 calculation", required=True)
args = vars(parser.parse_args())

gradvec = pQ.gradient(args["file"])
gradMAG = np.linalg.norm(gradvec, 'fro')

print('|force| = {:0.4f} eV bohr-1'.format(pQ.hartreetoeV(gradMAG)))
