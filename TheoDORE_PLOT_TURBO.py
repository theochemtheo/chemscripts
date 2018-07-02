#!/usr/bin/python3

'''
This is a template containing functions for reading the OmFrag.txt file and plotting the matrix
'''

import parseTURBO as pT
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt


# For reordering the matrices
def reorder(mat):
    reord = np.array(mat)
    reord = reord[:, reordlist]
    reord = reord[reordlist, :]
    return reord


# Use escf.out to get energies
escfFile = 'escf.out'
ExE = pT.cmtoeV(pT.excitations(escfFile))


# Extract the entire OmFragMat and keep only the lists of FragCT descriptors
OmFile = 'OmFrag.txt'
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
LLCT = np.ones([OmFragDim, OmFragDim], dtype=bool)
LLCT[0, 0] = False
LLCT = np.logical_xor(LLCT, bothMLCT)
LLCT = np.invert(np.logical_xor(LLCT, LCmask))

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


# Plotting
fig, ax = plt.subplots(5, 4, sharey='row', sharex='col')
for i in range(len(OmFragRaw)):
    OmFragMat = reorder(np.transpose(OmFragRaw[i].reshape([OmFragDim, OmFragDim])))
    stateE = ExE[i]
    # MC
    ax[np.floor_divide(i, 4), np.mod(i, 4)].pcolormesh(np.ma.array(OmFragMat, mask=reorder(MCmask)), cmap=plt.cm.Reds, shading='flat', edgecolors='black', vmin=0, vmax=1, linewidth=0.5)
    # MLCT/LMCT
    ax[np.floor_divide(i, 4), np.mod(i, 4)].pcolormesh(np.ma.array(OmFragMat, mask=reorder(bothMLCT)), cmap=plt.cm.Purples, shading='flat', edgecolors='black', vmin=0, vmax=1, linewidth=0.5)
    # LC
    ax[np.floor_divide(i, 4), np.mod(i, 4)].pcolormesh(np.ma.array(OmFragMat, mask=reorder(LCmask)), cmap=plt.cm.Greens, shading='flat', edgecolors='black', vmin=0, vmax=1, linewidth=0.5)
    # LLCT
    ax[np.floor_divide(i, 4), np.mod(i, 4)].pcolormesh(np.ma.array(OmFragMat, mask=reorder(LLCT)), cmap=plt.cm.Blues, shading='flat', edgecolors='black', vmin=0, vmax=1, linewidth=0.5)
    ax[np.floor_divide(i, 4), np.mod(i, 4)].set(aspect='equal', xticks=tickpos, yticks=tickpos)
    ax[np.floor_divide(i, 4), np.mod(i, 4)].tick_params(axis='both', which='both', bottom=False, left=False)
    ax[np.floor_divide(i, 4), np.mod(i, 4)].set_xticklabels(reordlab, rotation=90, horizontalalignment='center')
    ax[np.floor_divide(i, 4), np.mod(i, 4)].set_yticklabels(reordlab, verticalalignment='center')
    ax[np.floor_divide(i, 4), np.mod(i, 4)].set_title('{:1.2f} eV'.format(stateE), loc='center')
    for (x, y), z in np.ndenumerate(OmFragMat):
        if 0.7 > z >= 0.01:
            ax[np.floor_divide(i, 4), np.mod(i, 4)].text(y + 0.5, x + 0.5, '{:0.0f}'.format(z * 100), ha='center', va='center')
        elif z >= 0.7:
            ax[np.floor_divide(i, 4), np.mod(i, 4)].text(y + 0.5, x + 0.5, '{:0.0f}'.format(z * 100), ha='center', va='center', color='white')
    for each in range(OmFragDim):
        ax[np.floor_divide(i, 4), np.mod(i, 4)].text(each + 0.5, OmFragDim + 0.5, '{:0.0f}'.format(sum(OmFragMat[:, each]) * 100), ha='center', va='center')
        ax[np.floor_divide(i, 4), np.mod(i, 4)].text(OmFragDim + 0.5, each + 0.5, '{:0.0f}'.format(sum(OmFragMat[each, :]) * 100), ha='center', va='center')
plt.subplots_adjust(left=0.07, right=0.96, bottom=0.07, top=0.95, wspace=0.36, hspace=0.18)
# plt.tight_layout()
plt.savefig('OmFragPlot.pdf')
plt.show()
