#!/usr/bin/python3

# This script adapted (with permission) from the script "PlotBand", a script by J. Grant Hill.

import numpy as np
import matplotlib.pyplot as plt
import cclib
import argparse
import lineshapes


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
    thispeak = lineshapes.IRLorentzian(x, peak, parsedfile.vibirs[count], args["fwhm"])
    composite += thispeak

# If saving has been asked for
if args["save"] is not None:
    with open(args["save"], 'w') as csv:
        # Print a pre-amble so that csv is human readable
        csv.write("# Extracted & convoluted IR spectrum, with Lorentzian lineshapes, FWHM = {} cm-1\n".format(args["fwhm"]))
        csv.write("# x / cm-1; epsilon / L mol-1 cm-1; vib. freq. / cm-1; intensity / km mol-1\n")
        for each in range(max(len(x), len(parsedfile.vibfreqs))):
            try:
                csv.write("{},{},{},{}\n".format(x[each], composite[each], parsedfile.vibfreqs[each], parsedfile.vibirs[each]))
            except IndexError:
                if len(x) > len(parsedfile.vibfreqs):
                    csv.write("{},{}\n".format(x[each], composite[each]))
                else:
                    csv.write("{},{}\n".format(parsedfile.vibfreqs[each], parsedfile.vibirs[each]))

# If not in quiet mode
if args["quiet"] is None:
    # Set up the plot
    fig, ax1 = plt.subplots()
    # Axis 1 is the convoluted IR spectrum
    ax1.plot(x, composite)
    plt.xlabel('Frequency / cm$^{-1}$')
    ax1.set_ylabel('Molar absorption coefficient, $\epsilon$ / L mol$^{-1}$ cm$^{-1}$')
    # Axis 2 is the stick spectrum, which shares the x axis
    ax2 = ax1.twinx()
    ax2.vlines(parsedfile.vibfreqs, 0, parsedfile.vibirs)
    ax2.set_ylabel('IR intensity / km mol$^{-1}$')
    # Change formatting and plot
    fig.tight_layout()
    plt.show()
