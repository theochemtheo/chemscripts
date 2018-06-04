#!/usr/bin/python3

'''
This script uses cclib to extract vibrational normal mode vectors from an output file and prints data to the screen
'''


def VIBbilder(VIBMODE):
    '''
    This function returns the formatted .bild for the input vibrational mode

    reads in a vector from parseQCHEM

    returns tuple containing strings with terminated with newlines

    '''

    VIBMODEbild = ()

    # Vector settings
    arrow_color = "light sea green"
    # Scales magnitude of vector
    scale_fact = np.float(args["mult"])
    # radius of the stick. This should be >1
    arrow_stk_rad = ".1"
    # radius of the base of the tip cone. This should be ~4*arrow_stk_rad
    arrow_tip_rad = ".2"
    # fraction of total vector length taken up by stick
    arrow_stk_frc = "0.75"

    # Get the atom coordinates
    coords = infile.atomcoords[0]

    # Pre-amble
    VIBMODEbild = VIBMODEbild + (".color {}".format(arrow_color), )
    # For each atom, both forwards and backwards pointing vectors
    for atom in range(infile.natom):
        # Omit very short vectors
        if np.linalg.norm(VIBMODE[atom]) > args["cut"]:
            vec_start = coords[atom]
            vec_end_pos = vec_start + np.multiply(VIBMODE[atom], scale_fact)
            vec_end_neg = vec_start - np.multiply(VIBMODE[atom], scale_fact)
            VIBMODEbild = VIBMODEbild + (".arrow {:2.6f} {:2.6f} {:2.6f} {:2.6f} {:2.6f} {:2.6f} {} {} {}".format(vec_start[0], vec_start[1], vec_start[2], vec_end_pos[0], vec_end_pos[1], vec_end_pos[2], arrow_stk_rad, arrow_tip_rad, arrow_stk_frc), )
            VIBMODEbild = VIBMODEbild + (".arrow {:2.6f} {:2.6f} {:2.6f} {:2.6f} {:2.6f} {:2.6f} {} {} {}".format(vec_start[0], vec_start[1], vec_start[2], vec_end_neg[0], vec_end_neg[1], vec_end_neg[2], arrow_stk_rad, arrow_tip_rad, arrow_stk_frc), )

    return VIBMODEbild


import argparse
import numpy as np
import cclib
import os.path


parser = argparse.ArgumentParser(description="This script vibrational normal modes from a Gaussian 09 file and prints them to the screen in the Chimera .bild format")
parser.add_argument("-i", dest="file", metavar="file", help="input Gaussian 09 calculation", required=True)
parser.add_argument("-q", dest="quiet", help="don't print to screen", required=False, default=False, action='store_true')
parser.add_argument("-s", dest="save", help="save to file, with name <input file>.m<N>.bild, where N is the number of the mode", default=False, action='store_true')
parser.add_argument("-c", dest="cut", metavar="cut", help="threshold for atomic contributions, as a fraction. Default = 0.02. (very short vectors look bad)", default=0.02, required=False)
parser.add_argument("-m", dest="mult", metavar="mult", help="multiplicative factor for increasing the size of the vector. Default = 1.", default=1, required=False)
args = vars(parser.parse_args())

infile = cclib.parser.ccopen(args["file"]).parse()

# Determine padding requirements for output name
n_modes = len(infile.vibfreqs)
pad_length = len(str(n_modes))

for index, mode in enumerate(infile.vibdisps):
    bild = VIBbilder(mode)
    if not args["quiet"]:
        for line in range(len(bild)):
            print('{}'.format(bild[line]))
    if args["save"]:
        basename = os.path.basename(args["file"])[:-4]
        thisname = '{0}.m{1:0{1}d}.bild'.format(basename, index + 1, pad_length)
        if not os.path.exists(thisname):
            with open(thisname, 'a') as bild_file:
                for line in range(len(bild)):
                    bild_file.write('{}\n'.format(bild[line]))
        else:
            print('Warning, {} already exists!'.format(thisname))
