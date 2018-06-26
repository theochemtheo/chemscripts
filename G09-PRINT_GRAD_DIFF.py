#!/usr/bin/python3

'''
This script uses parseGAUSSIAN to extract gradient vectors from an output file and prints data to the screen
'''


import argparse
import cclib
import numpy as np
import parseGAUSSIAN as pG


parser = argparse.ArgumentParser(description="This script prints gradients from a Gaussian 09 file")
parser = argparse.ArgumentParser(description="This script extracts gradient differences from two G09 files and prints them to the screen in the Chimera .bild format")
parser.add_argument("-i", dest="state1", metavar="file 1", help="input G09 calculation containing gradient 1", required=True)
parser.add_argument("-j", dest="state2", metavar="file 2", help="input G09 calculation containing gradient 2", required=True)
parser.add_argument("-s", dest="show", help="print the vector to the screen", required=False, default=False, action='store_true')
args = vars(parser.parse_args())

infileI = cclib.parser.ccopen(args["state1"]).parse()
infileJ = cclib.parser.ccopen(args["state2"]).parse()

gradvecI = pG.gradient(args["state1"])
gradvecJ = pG.gradient(args["state2"])

graddiffvec = gradvecI - gradvecJ
graddiffMAG = np.linalg.norm(graddiffvec)

print('|force| = {:0.4f} eV bohr-1'.format(pG.hartreetoeV(graddiffMAG)))

if args["show"]:
    print('')
    print('{}'.format(graddiffvec))
