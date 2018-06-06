#!/usr/bin/python3

'''
This script uses parseQCHEM to extract gradient vectors from an output file and prints data to the screen
'''


import argparse
import numpy as np
import parseQCHEM as pQ


parser = argparse.ArgumentParser(description="This script prints gradients from a Q-Chem 5.0 file")
parser.add_argument("-i", dest="file", metavar="file", help="input Q-Chem 5.0 calculation", required=True)
parser.add_argument("-s", dest="show", help="print the vector to the screen", required=False, default=False, action='store_true')
args = vars(parser.parse_args())

gradvec = pQ.gradient(args["file"])
gradMAG = np.linalg.norm(gradvec, 'fro')

print('|force| = {:0.4f} eV bohr-1'.format(pQ.hartreetoeV(gradMAG)))

if args["show"]:
    print('')
    print('{}'.format(gradvec))
