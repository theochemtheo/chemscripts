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

NACVS = {}
DATA = {}
NACVS, DATA = pQ.NACV(args["file"])

for i in NACVS.keys():
    NACmag = np.linalg.norm(NACVS[i], 'fro')
    DCwithETF = NACmag / pQ.hartreetoeV(DATA['deltaE_{}'.format(i)])
    print('{}: delta E = {:0.3f} eV,  |NAC force| = {:0.4f} eV bohr-1, |Deriv. Coupling| = |NAC|/deltaEij = {:0.4f} bohr-1'.format(i, pQ.hartreetoeV(DATA['deltaE_{}'.format(i)]), NACmag, DCwithETF))
