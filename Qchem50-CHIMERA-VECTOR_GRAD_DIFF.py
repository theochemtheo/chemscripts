#!/usr/bin/python3

'''
This script uses parseQCHEM to extract non-adiabatic coupling vectors from an output file and prints data to the screen
'''


def GRADbilder(GRAD):
    '''
    This function returns the formatted .bild for the input NAC

    reads in a vector from parseQCHEM

    returns tuple containing strings with terminated with newlines

    '''

    # Norm the GRAD
    GRADmag = np.linalg.norm(GRAD, 'fro')
    nGRAD = np.divide(GRAD, GRADmag)

    GRADbild = ()

    # Vector settings
    arrow_color = "forest green"
    # Scales magnitude of vector
    scale_fact = args["mult"]
    # radius of the stick. This should be >1
    arrow_stk_rad = ".1"
    # radius of the base of the tip cone. This should be ~4*arrow_stk_rad
    arrow_tip_rad = ".2"
    # fraction of total vector length taken up by stick
    arrow_stk_frc = "0.75"

    # Get the atom coordinates
    coords = infileI.atomcoords[0]

    # Create the t
    GRADbild = GRADbild + (".color {}".format(arrow_color), )
    # For each atom, add a line to the tuple
    for atom in range(infileI.natom):
        # Omit very short vectors
        if np.linalg.norm(nGRAD[atom]) > args["cut"]:
            vec_start = coords[atom]
            vec_end = vec_start + np.multiply(nGRAD[atom], scale_fact)
            # GRADbild = GRADbild + (".arrow {} {} {} {} {} {}".format(vec_start[0], vec_start[1], vec_start[2], vec_end[0], vec_end[1], vec_end[2]), )
            GRADbild = GRADbild + (".arrow {:2.6f} {:2.6f} {:2.6f} {:2.6f} {:2.6f} {:2.6f} {} {} {}".format(vec_start[0], vec_start[1], vec_start[2], vec_end[0], vec_end[1], vec_end[2], arrow_stk_rad, arrow_tip_rad, arrow_stk_frc), )

    return GRADbild


import argparse
import numpy as np
import parseQCHEM as pQ
import cclib
import os.path


parser = argparse.ArgumentParser(description="This script extracts gradient differences from two Q-Chem 5.0 files and prints them to the screen in the Chimera .bild format")
parser.add_argument("-i", dest="state1", metavar="file 1", help="input Q-Chem 5.0 calculation containing gradient 1", required=True)
parser.add_argument("-j", dest="state2", metavar="file 2", help="input Q-Chem 5.0 calculation containing gradient 2", required=True)
parser.add_argument("-q", dest="quiet", help="don't print to screen", required=False, default=False, action='store_true')
parser.add_argument("-s", dest="name", help="save to this output file with suffix <name>.graddiff.bild", default=False)
parser.add_argument("-c", dest="cut", metavar="cut", help="threshold for atomic contributions, as a fraction. Default = 0.02. (very short vectors look bad)", default=0.02, required=False)
parser.add_argument("-m", dest="mult", metavar="mult", help="multiplicative factor for increasing the size of the vector. Default = 2.", default=2, required=False)
args = vars(parser.parse_args())

infileI = cclib.parser.ccopen(args["state1"]).parse()
infileJ = cclib.parser.ccopen(args["state2"]).parse()

gradvecI = pQ.gradient(args["state1"])
gradvecJ = pQ.gradient(args["state2"])

graddiffvec = gradvecI - gradvecJ

bild = GRADbilder(graddiffvec)
if not args["quiet"]:
    for line in range(len(bild)):
        print('{}'.format(bild[line]))
if args["name"]:
    basename = os.path.basename(args["name"])
    thisname = '{}.graddiff.bild'.format(basename)
    if not os.path.exists(thisname):
        with open(thisname, 'a') as bild_file:
            for line in range(len(bild)):
                bild_file.write('{}\n'.format(bild[line]))
    else:
        print('Warning, {} already exists!'.format(thisname))
