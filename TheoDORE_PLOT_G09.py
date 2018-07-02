#!/usr/bin/python3

'''
This is a template containing functions for reading the OmFrag.txt file and plotting the matrix
'''

import parseFUNCS as pF
import numpy as np
import argparse
import cclib
import matplotlib as mpl
import matplotlib.pyplot as plt


# For reordering the matrices
def reorder(mat):
    reord = np.array(mat)
    reord = reord[:, reordlist]
    reord = reord[reordlist, :]
    return reord


parser = argparse.ArgumentParser(description="This script generates nice electron-hole correlation plots")
parser.add_argument(dest="input", metavar="filename", help="input Gaussian 09 calculation to be analysed")
parser.add_argument("-o", dest="omfrag", metavar="filename", help="related OmFrag file from TheoDORE. If none is supplied, <input file>.OmFrag is assumed", required=False, default=False)
parser.add_argument("-q", dest="quiet", help="by default the correlation plot is shown to the screen. Use -q to prevent this", required=False, default=False, action='store_true')
parser.add_argument("-s", dest="save", help="save the correlation plot as \'<input file>_OmFragPlot.pdf\'", required=False, default=False, action='store_true')
parser.add_argument("-c", dest="col", metavar="columns", help="number of columns of states to plot (default=4)", required=False, default=4, type=int)
parser.add_argument("-r", dest="row", metavar="rows", help="number of rows of states to plot (default=5)", required=False, default=5, type=int)
args = vars(parser.parse_args())


# Basename for later
basename = args["input"][:-4]

# Use escf.out to get energies
parsed = cclib.parser.ccopen('{}'.format(args["input"])).parse()
ExE = pF.cmtoeV(parsed.etenergies)
ExS = parsed.etsyms


# Extract the entire OmFragMat and keep only the lists of FragCT descriptors
if args["omfrag"]:
    OmFile = '{}'.format(args["omfrag"])
else:
    OmFile = '{}.OmFrag'.format(basename)
OmFragRaw = np.genfromtxt(OmFile, skip_header=1)[:, 2:]
# Dimensionality of matrix
OmFragDim = int(np.sqrt(OmFragRaw.shape[1]))

# Create some masks for silly plotting
# MC
MCmask = np.ones([OmFragDim, OmFragDim], dtype=bool)
MCmask[0, 0] = False
# MLCT and LMCT
MLCTmask = np.ones([OmFragDim, OmFragDim], dtype=bool)
MLCTmask[0, 1:] = False
LMCTmask = np.ones([OmFragDim, OmFragDim], dtype=bool)
LMCTmask[1:, 0] = False
bothMLCT = np.invert(np.logical_xor(MLCTmask, LMCTmask))
# LC
LCmask = np.invert(np.eye(OmFragDim, dtype=bool))
LCmask[0, 0] = True
# LLCT
LLCTmask = np.ones([OmFragDim, OmFragDim], dtype=bool)
LLCTmask[0, 0] = False
LLCTmask = np.logical_xor(LLCTmask, bothMLCT)
LLCTmask = np.invert(np.logical_xor(LLCTmask, LCmask))

# Some sensible defaults
axlabels = ["M"]
for x in range(OmFragDim - 1):
    label = 'L{}'.format(x + 1)
    axlabels.append(label)

tickpos = [x + 0.5 for x in range(OmFragDim)]
reordlist = [x for x in range(OmFragDim)]

# Place custom settings for labels etc. here!


# Reorder the labels
reordlab = [axlabels[i] for i in reordlist]

# Set up mpl
mpl.rcParams['figure.figsize'] = [7, 8]
mpl.rcParams['figure.dpi'] = 100
mpl.rcParams['savefig.dpi'] = 300
mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.unicode'] = True
mpl.rcParams['text.latex.preamble'] = [
    r'\usepackage{amsmath,amssymb,amsfonts,textcomp,array}',
    r'\usepackage{txfonts}',
    r'\usepackage{sansmath}',
    r'\renewcommand*\sfdefault{phv}',
    r'\renewcommand{\familydefault}{\sfdefault}',
    r'\sansmath']

# Change border box linewidth
mpl.rcParams['axes.linewidth'] = 1

# plottings variables
plot_cols = args["col"]
plot_rows = args["row"]
n_states = len(ExE)
n_plot = plot_rows * plot_cols
if n_plot > n_states:
    plot_rows = np.int(np.floor(n_states / plot_cols)) + 1

# Plotting
fig, ax = plt.subplots(plot_rows, plot_cols, sharey='row', sharex='col')
for i in range(n_states):
    ThisOmFragMat = reorder(np.transpose(OmFragRaw[i].reshape([OmFragDim, OmFragDim])))
    stateE = ExE[i]
    # MC
    ax[np.floor_divide(i, plot_cols), np.mod(i, plot_cols)].pcolormesh(np.ma.array(ThisOmFragMat, mask=reorder(MCmask)), cmap=plt.cm.Reds, shading='flat', edgecolors='black', vmin=0, vmax=1, linewidth=0.5)
    # MLCT/LMCT
    ax[np.floor_divide(i, plot_cols), np.mod(i, plot_cols)].pcolormesh(np.ma.array(ThisOmFragMat, mask=reorder(bothMLCT)), cmap=plt.cm.Purples, shading='flat', edgecolors='black', vmin=0, vmax=1, linewidth=0.5)
    # LC
    ax[np.floor_divide(i, plot_cols), np.mod(i, plot_cols)].pcolormesh(np.ma.array(ThisOmFragMat, mask=reorder(LCmask)), cmap=plt.cm.Greens, shading='flat', edgecolors='black', vmin=0, vmax=1, linewidth=0.5)
    # LLCT
    ax[np.floor_divide(i, plot_cols), np.mod(i, plot_cols)].pcolormesh(np.ma.array(ThisOmFragMat, mask=reorder(LLCTmask)), cmap=plt.cm.Blues, shading='flat', edgecolors='black', vmin=0, vmax=1, linewidth=0.5)
    ax[np.floor_divide(i, plot_cols), np.mod(i, plot_cols)].set(aspect='equal', xticks=tickpos, yticks=tickpos)
    ax[np.floor_divide(i, plot_cols), np.mod(i, plot_cols)].tick_params(axis='both', which='both', bottom=False, left=False)
    ax[np.floor_divide(i, plot_cols), np.mod(i, plot_cols)].set_xticklabels(reordlab, rotation=90, horizontalalignment='center')
    ax[np.floor_divide(i, plot_cols), np.mod(i, plot_cols)].set_yticklabels(reordlab, verticalalignment='center')
    ax[np.floor_divide(i, plot_cols), np.mod(i, plot_cols)].set_title('{}, {:1.2f} eV'.format(ExS[i][0:7], stateE), loc='center', y=1.1)
    for (x, y), z in np.ndenumerate(ThisOmFragMat):
        if 0.7 > z >= 0.01:
            ax[np.floor_divide(i, plot_cols), np.mod(i, plot_cols)].text(y + 0.5, x + 0.5, '{:0.0f}'.format(z * 100), ha='center', va='center')
        elif z >= 0.7:
            ax[np.floor_divide(i, plot_cols), np.mod(i, plot_cols)].text(y + 0.5, x + 0.5, '{:0.0f}'.format(z * 100), ha='center', va='center', color='white')
    for each in range(OmFragDim):
        ax[np.floor_divide(i, plot_cols), np.mod(i, plot_cols)].text(each + 0.5, OmFragDim + 0.5, '{:0.0f}'.format(sum(ThisOmFragMat[:, each]) * 100), ha='center', va='center')
        ax[np.floor_divide(i, plot_cols), np.mod(i, plot_cols)].text(OmFragDim + 0.5, each + 0.5, '{:0.0f}'.format(sum(ThisOmFragMat[each, :]) * 100), ha='center', va='center')
plt.subplots_adjust(left=0.07, right=0.96, bottom=0.07, top=0.95, wspace=0.36, hspace=0.18)
# plt.tight_layout()
if args["save"]:
    plt.savefig('{}_OmFragPlot.pdf'.format(basename))
if not args["quiet"]:
    plt.show()
