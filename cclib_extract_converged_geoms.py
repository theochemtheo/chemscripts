#!/usr/bin/env python3

"""
This script is based on cclib_extract_geom.py, written by Eric Berquist and published on https://chemistry.stackexchange.com/a/62562
"""


from __future__ import print_function
import os.path
from cclib.io import ccread, ccwrite
import logging
import numpy as np


def getargs():
    """Get command-line arguments."""

    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('outputfilename', nargs='+')

    parser.add_argument('--trajectory',
                        action='store_true',
                        help="""Extract all possible geometries into a trajectory-like file?""")
    parser.add_argument('--suffix')

    args = parser.parse_args()

    return args


if __name__ == '__main__':

    args = getargs()

    for outputfilename in args.outputfilename:

        data = ccread(outputfilename, loglevel=logging.ERROR)
        converged_geoms = np.where(data.optstatus == 4)[0]

        stub = os.path.splitext(outputfilename)[0]
        if args.suffix:
            xyzfilename_prefix = '{}.{}'.format([stub, args.suffix])
        else:
            xyzfilename_prefix = stub

        for step, index in enumerate(converged_geoms):
            xyzfilename = '{}_{}.xyz'.format(xyzfilename_prefix, step + 1)
            ccwrite(data, outputdest=xyzfilename, indices=index, outputtype='xyz')

