#!/usr/bin/python3

# This script adapted (with permission) from the script "PlotBand", a script by J. Grant Hill.

import sys
# Check for the required modules, try to exit gracefully if not found
import imp
try:
    imp.find_module('numpy')
    foundnp = True
except ImportError:
    foundnp = False
try:
    imp.find_module('matplotlib')
    foundplot = True
except ImportError:
    foundplot = False
try:
    imp.find_module('cclib')
    foundcclib = True
except ImportError:
    foundcclib = False
try:
    imp.find_module('argparse')
    foundarg = True
except ImportError:
    foundarg = False
if not foundnp:
    print("Numpy is required. Exiting")
    sys.exit()
if not foundplot:
    print("Matplotlib is required. Exiting")
    sys.exit()
if not foundcclib:
    print("cclib is required. Exiting")
    sys.exit()
if not foundarg:
    print("Argh! Argparse is required. Exiting")
    sys.exit()
import numpy as np
import matplotlib.pyplot as plt
import cclib
import argparse


def Lorentzian(x, band, strength, FWHM):
    # Lorentzian lineshape function
    # This is based on a document by Jens Spanget-Larsen, which may be found at dx.doi.org/10.13140/RG.2.1.4181.6160
    HWHM = FWHM / 2
    bandshape = ((100 / np.log(10) / np.pi * strength) / HWHM) / (1 + ((x - band) / HWHM)**2)
    return bandshape


# Set up argument parsing
parser = argparse.ArgumentParser(description="This script uses cclib to extract a Lorentzian broadened IR spectrum from a Gaussian 09 calculation")
parser.add_argument("-i", dest="file", metavar="file", help="input Gaussian 09 frequency calculation", required=True)
parser.add_argument("-f", dest="fwhm", metavar="fwhm", help="the FWHM in cm-1 of the desired Lorentzians. Default: 6 cm-1", type=int, required=False, default=6)
parser.add_argument("-b", dest="begin", metavar="begin", help="beginning (in cm-1) of plot. Default: 0 cm-1", type=int, required=False, default=0)
parser.add_argument("-e", dest="end", metavar="end", help="end (in cm-1) of plot. Default: 4000 cm-1", type=int, required=False, default=4000)
parser.add_argument("-p", dest="points", metavar="points", help="number of points in plot. Default: 8000", type=int, required=False, default=8000)
parser.add_argument("-s", dest="save", metavar="save", help="file in which output will be save (optional)", required=False, default=None)
parser.add_argument("-q", dest="quiet", metavar="quiet", help="by default a matplotlib window will appear, use -q True to prevent this", required=False, default=None)
args = vars(parser.parse_args())

# Parse the file using cclib
parsedfile = cclib.io.ccread(args["file"])

# Create the x axis using the settings
x = np.linspace(args["begin"], args["end"], args["points"])

# Make the Lorentzians and add them together
composite = 0
for count, peak in enumerate(parsedfile.vibfreqs):
    thispeak = Lorentzian(x, peak, parsedfile.vibirs[count], args["fwhm"])
    composite += thispeak

# If
if args["quiet"] is None:
    # Set up the plot
    fig, ax1 = plt.subplots()
    ax1.plot(x, composite)
    plt.xlabel('Frequency / cm$^{-1}$')
    ax1.set_ylabel('Molar absorption coefficient, $\epsilon$ / L mol$^{-1}$ cm$^{-1}$')
    ax2 = ax1.twinx()
    ax2.vlines(parsedfile.vibfreqs, 0, parsedfile.vibirs)
    ax2.set_ylabel('IR intensity / km mol$^{-1}$')
    fig.tight_layout()
    plt.show()

if args["save"] is not None:
    print("blargh")
