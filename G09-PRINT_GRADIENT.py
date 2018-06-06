#!/usr/bin/python3

'''
This script uses parseGAUSSIAN to extract gradient vectors from an output file and prints data to the screen
'''


import argparse
import numpy as np
import parseGAUSSIAN as pG


parser = argparse.ArgumentParser(description="This script prints gradients from a Gaussian 09 file")
parser.add_argument("-i", dest="file", metavar="file", help="input Gaussian 09 calculation", required=True)
parser.add_argument("-s", dest="show", help="print the vector to the screen", required=False, default=False, action='store_true')
args = vars(parser.parse_args())

gradvec = pG.gradient(args["file"])
gradMAG = np.linalg.norm(gradvec, 'fro')

print('|force| = {:0.4f} eV bohr-1'.format(pG.hartreetoeV(gradMAG)))

if args["show"]:
    print('')
    print('{}'.format(gradvec))
