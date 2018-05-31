#!/usr/bin/python3

'''
This script uses parseQCHEM to extract non-adiabatic coupling vectors from an output file and prints data to the screen
'''


import argparse
import numpy as np
import parseQCHEM as pQ


parser = argparse.ArgumentParser(description="This script prints gradient difference vectors from two Q-Chem 5.0 files")
parser.add_argument("-i", dest="state1", metavar="file 1", help="input Q-Chem 5.0 calculation containing gradient 1", required=True)
parser.add_argument("-j", dest="state2", metavar="file 2", help="input Q-Chem 5.0 calculation containing gradient 2", required=True)
args = vars(parser.parse_args())

gradvecI = pQ.gradient(args["state1"])
gradvecJ = pQ.gradient(args["state2"])
gradMAGI = np.linalg.norm(gradvecI, 'fro')
gradMAGJ = np.linalg.norm(gradvecJ, 'fro')

graddiffvec = gradvecI - gradvecJ
graddiffMAG = np.linalg.norm(graddiffvec, 'fro')

print('|force i| = {:0.4f} eV bohr-1, |force j| = {:0.4f} eV bohr-1, |force i-j| = {:0.4f} eV bohr-1'.format(pQ.hartreetoeV(gradMAGI), pQ.hartreetoeV(gradMAGJ), pQ.hartreetoeV(graddiffMAG)))

print('gradient difference vector in hartree bohr-1:')
print('{}'.format(graddiffvec))
