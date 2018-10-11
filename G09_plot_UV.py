#!/usr/bin/python3

# This script adapted (with permission) from the script "PlotBand", a script by J. Grant Hill.

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import cclib
import argparse


# def UVLorentzian(x, band, strength, FWHM):
#     # Lorentzian lineshape function with UV prefactor
#     # This is based on a document by Jens Spanget-Larsen, which may be found at dx.doi.org/10.13140/RG.2.1.4181.6160
#     HWHM = FWHM / 2
#     bandshape = (2.175e8) * (strength / FWHM) / (1 + ((x - band) / HWHM)**2)
#     return bandshape


def UVGaussian(x, band, strength, FWHM):
    # Gaussian lineshape function with UV prefactor
    bandshape = (2.175e8) * (strength / FWHM) * np.exp(-2.772 * ((x - band) / FWHM)**2)
    return bandshape


def tenthousands(x, pos):
    return '%1d' % (x * 1e-3)


# Set up argument parsing
parser = argparse.ArgumentParser(description="This script uses cclib to plot a broadened UV/vis spectrum from a Gaussian 09 calculations")
parser.add_argument("-i", dest="file", metavar="file", help="input Gaussian 09 calculation", required=True)
parser.add_argument("-f", dest="fwhm", metavar="FWHM", help="the FWHM in cm-1 of the desired Gaussians. Default: 1500 cm-1", type=int, required=False, default=1500)
parser.add_argument("-b", dest="begin", metavar="begin", help="beginning (in nm) of plot. Default: 300 nm", type=int, required=False, default=300)
parser.add_argument("-e", dest="end", metavar="end", help="end (in nm) of plot. Default: 800 nm", type=int, required=False, default=800)
parser.add_argument("-p", dest="points", metavar="points", help="number of points in plot. Default: 8000", type=int, required=False, default=8000)
parser.add_argument("-s", dest="save", metavar="save", help="file in which plot will be saved (optional)", required=False, default=None)
parser.add_argument("-x", dest="shift", metavar="shift", help="the transition energies can be shifted by an amount in cm-1. Default: 0", required=False, default=0.0, type=float)
parser.add_argument("-d", dest="data", help="should the sticks and plot data be saved? (optional)", required=False, default=None, action='store_true')
parser.add_argument("-q", dest="quiet", help="by default a matplotlib window will appear, use -q to prevent this", required=False, default=None, action='store_true')
parser.add_argument("-a", "--no-latex", dest="allergy", help="disable plotting with LaTeX", required=False, default=None, action='store_true')
args = vars(parser.parse_args())

# Create the x axis using the settings
x = np.linspace(1e7 / args["begin"], 1e7 / args["end"], args["points"])
smoothx = np.linspace(1e7 / args["begin"], 1e7 / args["end"], 150)

# Files to be parsed
parsed_input = cclib.parser.ccopen(args["file"]).parse()

# If a shift has been applied, use it
parsed_input.etenergies = parsed_input.etenergies + args["shift"]

# Basename
base_fname = args["file"].replace('.log', '')

# Peak fitting: make a np array of the correct dimension first
composite_G_abs = np.zeros([1, args["points"]])
for count, peak in enumerate(parsed_input.etenergies):
    thispeak = UVGaussian(x, peak, parsed_input.etoscs[count], args["fwhm"])
    composite_G_abs += thispeak
# Have to transpose
composite_G_abs = np.transpose(composite_G_abs)

if args["data"] is True:
    # Include shift in header if necessary
    if args["shift"] != 0.0:
        stk_header = 'Excitation energies / cm-1 (shifted by {} cm-1) and oscillator strengths / dimensionless, extracted from {}'.format(args["shift"], args["file"])
        plt_header = 'Convoluted UV/vis spectrum with Gaussians, FWHM = {}, (transitions shifted by {} cm-1) from file {}.\nx = wavelength / nm\t\ty = epsilon / L mol-1 cm-1'.format(args["fwhm"], args["shift"], args["file"])
    else:
        stk_header = 'Excitation energies / cm-1 and oscillator strengths / dimensionless, extracted from {}'.format(args["file"])
        plt_header = 'Convoluted UV/vis spectrum with Gaussians, FWHM = {}, from file {}.\nx = wavelength / nm\t\ty = epsilon / L mol-1 cm-1'.format(args["fwhm"], args["file"])
    stk_array = np.column_stack((parsed_input.etenergies, parsed_input.etoscs))
    plt_array = np.column_stack((1e7 / x, composite_G_abs))
    np.savetxt('{}.stk'.format(base_fname), stk_array, delimiter='\t', header=stk_header)
    np.savetxt('{}.plt'.format(base_fname), plt_array, delimiter='\t', header=plt_header)

# Set up mpl
mpl.rcParams['figure.figsize'] = [9, 6]
mpl.rcParams['figure.dpi'] = 100
mpl.rcParams['savefig.dpi'] = 300
if args["allergy"] is None:
    mpl.rcParams['text.usetex'] = True
    mpl.rcParams['text.latex.unicode'] = True
    mpl.rcParams['text.latex.preamble'] = [
        r'\usepackage{amsmath,amssymb,amsfonts,textcomp,array}',
        r'\usepackage{txfonts}',
        r'\usepackage{sansmath}',
        r'\renewcommand*\sfdefault{phv}',
        r'\renewcommand{\familydefault}{\sfdefault}',
        r'\sansmath']
# For altering y axis
tenks = mpl.ticker.FuncFormatter(tenthousands)

# Prepare the axes
fig, ax = plt.subplots()
# Create the twinned axes for oscillator strength
axes = [ax, ax.twinx()]
# labels - check for LaTeX allergies
if args["allergy"] is None:
    axes[0].set_xlabel('Wavelength, $\lambdaup$ / nm', fontsize=14)
    axes[0].set_ylabel('Molar absorption coefficient, $\\varepsilonup$ / $10^3\\times$ L mol$^{-1}$ cm$^{-1}$', fontsize=14)
    axes[1].set_ylabel('Oscillator strength, $\\mathit{f}$ / dimensionless', fontsize=14)
else:
    axes[0].set_xlabel('Wavelength, lambda / nm', fontsize=14)
    axes[0].set_ylabel('Molar absorption coefficient, epsilon / 1000 * L mol-1 cm-1', fontsize=14)
    axes[1].set_ylabel('Oscillator strength, f / dimensionless', fontsize=14)
# axes[0] is the convoluted UV spectrum
axes[0].set_ylim([0, 1e5])
axes[0].yaxis.set_major_formatter(tenks)
axes[0].tick_params(labelsize=12)
# Calculated
axes[0].plot(1e7 / (x), composite_G_abs, color='r', alpha=0.5)
# axes[1] is the stick spectrum
axes[1].set_ylim([0, 1])
axes[1].tick_params(labelsize=12)
axes[1].set_xlim([args["begin"], args["end"]])
# Plot sticks using scaled frequencies
axes[1].vlines(1e7 / parsed_input.etenergies, 0, parsed_input.etoscs, color='r', alpha=0.5)
# Change formatting and plot
fig.tight_layout()
# If not in quiet mode
if args["quiet"] is None:
    plt.show()
if args["save"] is not None:
    plt.savefig(args["save"])
