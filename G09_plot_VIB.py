#!/usr/bin/python3

# This script adapted (with permission) from the script "PlotBand", a script by J. Grant Hill.

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import cclib
import argparse


_daltontokg = 1.660539040e-27
_kBoltzmann = 1.38064852e-23
_hPlanck = 6.62607004e-34
_cSpeedOfLight = 299792458


def _Lorentzian(x, band, FWHM):
    HWHM = FWHM / 2
    probability = 1 / (1 + ((x - band) / HWHM)**2)
    return probability


def _Gaussian(x, band, FWHM):
    probability = np.exp(-2.772 * ((x - band) / FWHM)**2)
    return probability


def _IRprefac(strength, FWHM):
    # This is based on a document by Jens Spanget-Larsen, which may be found at dx.doi.org/10.13140/RG.2.1.4181.6160
    HWHM = FWHM / 2
    prefac = ((100 / np.log(10) / np.pi * strength) / HWHM)
    return prefac


def _Ramanprefac(band, activity, excitation, temperature):
    # This is based on eqn. 1 of V. Krishnakumar et al., J. Mol. Struct. 2004, 702, 9
    # dx.doi.org/10.1016/j.molstruc.2004.06.004
    prefac_num = (100 * (1e7 / excitation - band))**4 * (1e-40 * activity / _daltontokg)
    prefac_denom = (100 * band) * (1 - np.exp(-(_hPlanck * _cSpeedOfLight * 100 * band) / (_kBoltzmann * temperature)))
    prefac = prefac_num / prefac_denom
    return prefac


def IRLorentzian(x, band, strength, FWHM):
    # Lorentzian lineshape function with IR prefactor
    epsilon = _IRprefac(strength, FWHM) * _Lorentzian(x, band, FWHM)
    return epsilon


def IRGaussian(x, band, strength, FWHM):
    # Gaussian lineshape function with IR prefactor
    epsilon = _IRprefac(strength, FWHM) * _Gaussian(x, band, FWHM)
    return epsilon


def RamanLorentzian(x, band, activity, excitation, temperature, FWHM):
    # Lorentzian lineshape function with Raman prefactor
    intensity = _Ramanprefac(band, activity, excitation, temperature) * _Lorentzian(x, band, FWHM)
    return intensity


def RamanGaussian(x, band, activity, excitation, temperature, FWHM):
    # Gaussian lineshape function with Raman prefactor
    intensity = _Ramanprefac(band, activity, excitation, temperature) * _Gaussian(x, band, FWHM)
    return intensity


def tenthousands(x, pos):
    return '%1d' % (x * 1e-3)


# Set up argument parsing
parser = argparse.ArgumentParser(description="This script uses cclib to plot a broadened UV/vis spectrum from a Gaussian 09 calculations")
parser.add_argument("-i", dest="file", metavar="file", help="input Gaussian 09 calculation", required=True)
parser.add_argument("-f", dest="fwhm", metavar="FWHM", help="the FWHM in cm-1 of the desired broadening functions. Default: 8 cm-1", type=int, required=False, default=8)
parser.add_argument("-b", dest="begin", metavar="begin", help="beginning (in cm-1) of plot. Default: 0 cm-1", type=int, required=False, default=0)
parser.add_argument("-e", dest="end", metavar="end", help="end (in cm-1) of plot. Default: 3500 cm-1", type=int, required=False, default=3500)
parser.add_argument("-p", dest="points", metavar="points", help="number of points in plot. Default: 10000", type=int, required=False, default=10000)
parser.add_argument("-s", dest="save", metavar="save", help="file in which plot will be saved (optional)", required=False, default=None)
parser.add_argument("-d", dest="data", help="should the sticks and plot data be saved? (optional)", required=False, default=None, action='store_true')
parser.add_argument("-q", dest="quiet", help="by default a matplotlib window will appear, use -q to prevent this", required=False, default=None, action='store_true')
parser.add_argument("-L", dest="Lorentzian", help="by default, a Gaussian broadening function is applied, use -L to change this to Lorentzian", required=False, default=False, action='store_true')
parser.add_argument("-R", dest="Raman", help="by default, the IR spectrum is plotted, use -R to change this to non-resonant Raman", required=False, default=False, action='store_true')
parser.add_argument("-T", dest="temperature", help="T used in Raman intensity calculations. Default: 298.15 K", required=False, default=298.15)
parser.add_argument("-X", dest="excitation", help="Excitation wavelength (in nm) used in non-resonant Raman. Default: 1064 nm (Nd:YAG)", required=False, default=1064.0, type=float)
parser.add_argument("-n", dest="normalize", help="Set this to normalize the y-axis to 1", required=False, default=False, action='store_true')
parser.add_argument("-a", "--no-latex", dest="allergy", help="disable plotting with LaTeX", required=False, default=None, action='store_true')
args = vars(parser.parse_args())

# Create the x axis using the settings
x = np.linspace(args["begin"], args["end"], args["points"])

# Files to be parsed
parsed_input = cclib.parser.ccopen(args["file"]).parse()

# Basename
base_fname = args["file"].replace('.log', '')

# Make sure correct functions will be used
if args["Raman"]:
    spectrum_type = 'non-resonant Raman Scattering'
    intensities = parsed_input.vibramans
    intensities_name = 'Raman activity / Angstrom^4 Dalton^-1'
    if args["normalize"]:
        signal_name = 'Normalized Raman intensity / arb. units'
    else:
        signal_name = 'Raman intensity / m kg^-1'
    if args["Lorentzian"]:
        broadening_function = RamanLorentzian
        broadening_function_name = 'Lorentzian'
    else:
        broadening_function = RamanGaussian
        broadening_function_name = 'Gaussian'
else:
    spectrum_type = 'Infrared absorption'
    intensities = parsed_input.vibirs
    intensities_name = 'IR intensity / km mol^-1'
    if args["normalize"]:
        signal_name = 'Normalized molar absorption coefficient (epsilon) / arb. units'
    else:
        signal_name = 'Molar absorption coefficient (epsilon) / L mol^-1 cm^-1'
    if args["Lorentzian"]:
        broadening_function = IRLorentzian
        broadening_function_name = 'Lorentzian'
    else:
        broadening_function = IRGaussian
        broadening_function_name = 'Gaussian'

# Peak fitting: make a np array of the correct dimension first
composite_spectrum = np.zeros([1, args["points"]])
for count, peak in enumerate(parsed_input.vibfreqs):
    if args["Raman"]:
        thispeak = broadening_function(x, peak, intensities[count], args["excitation"], args["temperature"], args["fwhm"])
    else:
        thispeak = broadening_function(x, peak, intensities[count], args["fwhm"])
    composite_spectrum += thispeak
# Have to transpose
composite_spectrum = np.transpose(composite_spectrum)

# Normalize if necessary
if args["normalize"]:
    composite_spectrum = composite_spectrum / np.max(composite_spectrum)

if args["data"] is True:
    stk_header = 'Vibrational energies / cm^-1 and {}, extracted from {}'.format(intensities_name, args["file"])
    plt_header = 'Convoluted {} spectrum with convoluted with {}s, FWHM = {}, from file {}.\nx = wavenumber / cm^-1\t\ty = {}'.format(spectrum_type, broadening_function_name, args["fwhm"], args["file"], signal_name)
    stk_array = np.column_stack((parsed_input.vibfreqs, intensities))
    plt_array = np.column_stack((1e7 / x, composite_spectrum))
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

# Prepare the axes
fig, ax = plt.subplots()
# Create the twinned axes for transition probability
axes = [ax, ax.twinx()]
# labels - check for LaTeX allergies
if args["allergy"] is None:
    axes[0].set_xlabel('Frequency, $\\tilde{\\nu}$ / cm$^{-1}$', fontsize=14)
    if args["Raman"]:
        if args["normalize"]:
            axes[0].set_ylabel('{}'.format(signal_name), fontsize=14)
        else:
            axes[0].set_ylabel('Raman intensity, $I$ / m kg$^{-1}$', fontsize=14)
        axes[1].set_ylabel('Raman activity / \\AA$^4$ Da$^{-1}$', fontsize=14)
    else:
        if args["normalize"]:
            axes[0].set_ylabel('{}'.format(signal_name), fontsize=14)
        else:
            axes[0].set_ylabel('Molar absorption coefficient, $\\varepsilonup$ / L mol$^{-1}$ cm$^{-1}$', fontsize=14)
        axes[1].set_ylabel('IR intensity / km mol$^{-1}$', fontsize=14)
else:
    axes[0].set_xlabel('Frequency, nu-tilde / cm^-1', fontsize=14)
    axes[0].set_ylabel('{}'.format(signal_name), fontsize=14)
    axes[1].set_ylabel('{}'.format(intensities_name), fontsize=14)
# axes[0] is the convoluted spectrum
axes[0].set_ylim([0, 1.2 * np.max(composite_spectrum)])
axes[0].tick_params(labelsize=12)
# Calculated
axes[0].plot(x, composite_spectrum, color='r', alpha=0.5)
# axes[1] is the stick spectrum
axes[1].set_ylim([0, 1.5 * np.max(intensities)])
axes[1].tick_params(labelsize=12)
axes[1].set_xlim([args["begin"], args["end"]])
# Plot sticks using scaled frequencies
axes[1].vlines(parsed_input.vibfreqs, 0, intensities, color='r', alpha=0.5)
# Change formatting and plot
fig.tight_layout()
# If not in quiet mode
if args["quiet"] is None:
    plt.show()
if args["save"] is not None:
    plt.savefig(args["save"])
