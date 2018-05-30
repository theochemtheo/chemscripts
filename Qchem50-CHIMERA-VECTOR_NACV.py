#!/usr/bin/python3

'''
This script uses parseQCHEM to extract non-adiabatic coupling vectors from an output file and prints data to the screen
'''


def NACVbilder(NACV):
    '''
    This function returns the formatted .bild for the input NAC

    reads in a vector from parseQCHEM

    returns tuple containing strings with terminated with newlines

    '''

    # Norm the NACV
    NACmag = np.linalg.norm(NACV, 'fro')
    nNACV = np.divide(NACV, NACmag)

    NACVbild = ()

    # Vector settings
    arrow_color = "hot pink"
    # Scales magnitude of vector
    scale_fact = args["mult"]
    # radius of the stick. This should be >1
    arrow_stk_rad = ".1"
    # radius of the base of the tip cone. This should be ~4*arrow_stk_rad
    arrow_tip_rad = ".2"
    # fraction of total vector length taken up by stick
    arrow_stk_frc = "0.75"

    # Get the atom coordinates
    coords = infile.atomcoords[0]

    # Create the t
    NACVbild = NACVbild + (".color {}".format(arrow_color), )
    # For each atom, add a line to the tuple
    for atom in range(infile.natom):
        # Omit very short vectors
        if np.linalg.norm(nNACV[atom]) > args["cut"]:
            vec_start = coords[atom]
            vec_end = vec_start + np.multiply(nNACV[atom], scale_fact)
            # NACVbild = NACVbild + (".arrow {} {} {} {} {} {}".format(vec_start[0], vec_start[1], vec_start[2], vec_end[0], vec_end[1], vec_end[2]), )
            NACVbild = NACVbild + (".arrow {:2.6f} {:2.6f} {:2.6f} {:2.6f} {:2.6f} {:2.6f} {} {} {}".format(vec_start[0], vec_start[1], vec_start[2], vec_end[0], vec_end[1], vec_end[2], arrow_stk_rad, arrow_tip_rad, arrow_stk_frc), )

    return NACVbild


import argparse
import numpy as np
import parseQCHEM as pQ
import cclib
import os.path


parser = argparse.ArgumentParser(description="This script extracts non-adiabatic couplings from a Q-Chem 5.0 file and prints them to the screen in the Chimera .bild format")
parser.add_argument("-i", dest="file", metavar="file", help="input Q-Chem 5.0 calculation", required=True)
parser.add_argument("-q", dest="quiet", help="don't print to screen", required=False, default=False, action='store_true')
parser.add_argument("-s", dest="save", help="save to file, with name <input file>.<MI>-to-<MJ>.bild, where M = multiplicity, I, J = state number", default=False, action='store_true')
parser.add_argument("-c", dest="cut", metavar="cut", help="threshold for atomic contributions, as a fraction. Default = 0.02. (very short vectors look bad)", default=0.02, required=False)
parser.add_argument("-m", dest="mult", metavar="mult", help="multiplicative factor for increasing the size of the vector. Default = 2.", default=2, required=False)
args = vars(parser.parse_args())

infile = cclib.parser.ccopen(args["file"]).parse()

NACVS = {}
DATA = {}
NACVS, DATA = pQ.NACV(args["file"])

for pair in NACVS.keys():
    bild = NACVbilder(NACVS[pair])
    if not args["quiet"]:
        for line in range(len(bild)):
            print('{}'.format(bild[line]))
    if args["save"]:
        basename = os.path.basename(args["file"])[:-4]
        thisname = '{}.{}.bild'.format(basename, pair)
        if not os.path.exists(thisname):
            with open(thisname, 'a') as bild_file:
                for line in range(len(bild)):
                    bild_file.write('{}\n'.format(bild[line]))
        else:
            print('Warning, {} already exists!'.format(thisname))
