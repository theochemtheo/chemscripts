from numpy import loadtxt, nan, array, append, float
from re import findall, compile
from parseFUNCS import *


'''

parseTURBO provides similar functionality to cclib, but for Turbomole calculations

It has been tested with Turbomole version 7.2.

'''


def atomcoords(file):
    '''
    Extract atomic coordinates from a Turbomole calculation

    returns a numpy array of atomic coordinates in angstrom
    '''
    # Set up constants and arrays
    autoangstr = 0.529177
    rawmatched = []
    finmatched = array([]).reshape(0, 3)
    found = False
    # Read file and keep only lines that match coordinate lines
    readfile = open(file)
    for line in readfile:
        if found is True:
            if not line.strip() == '':
                rawmatched = append(rawmatched, line)
            else:
                break
        if "atomic coordinates" in line:
            found = True
    readfile.close()
    # Clear up the strings from these lines and create a numpy array analagous to the one from cclib
    for line in rawmatched:
        coordline = array(line.split()[0:3]).astype(float)
        finmatched = append(finmatched, [coordline * autoangstr], axis=0)
    return finmatched


def gsenergy(file):
    '''
    Extract ground state energy in hartrees

    returns a string
    '''
    gs = loadtxt(findall('(?<=Total energy\:).*', open(file).read()))[0]
    return gs


def esenergies(file):
    '''
    Extract excited state absolute energies in hartrees. This only works for escf calculatons.

    returns an array of strings
    '''
    es = loadtxt(findall('(?<=Total energy\:).*', open(file).read()))[1:]
    return es


def excitations(file):
    '''
    Extract excitation energies in cm-1. This only works for escf calculations.

    returns an array of strings
    '''
    excitations = loadtxt(findall('(?<=Excitation energy \/ cm\^\(-1\)\:).*', open(file).read()))
    return excitations


def twoCcheck(file):
    '''
    Checks if calculation is 2C or not.

    returns boolean
    '''
    istwoC = False
    with open(file, 'r') as f:
        for line in f.readlines():
            if 'Two-component modus switched on ! ' in line:
                istwoC = True
    return istwoC


def _oneCoscvelrep(file):
    '''
    Extract oscillator strengths in velocity representation from non-2C escf calculation.

    returns array of strings
    '''
    velrep = loadtxt(findall('(?<=velocity representation\:).*', open(file).read()))[0::2]
    return velrep


def _oneCosclenrep(file):
    '''
    Extract oscillator strengths in length representation from non-2C escf calculation.

    returns array of strings
    '''
    lenrep = loadtxt(findall('(?<=length representation\:).*', open(file).read()))[0::2]
    return lenrep


def _oneCoscmixrep(file):
    '''
    Extract oscillator strengths in mixed representation from non-2C escf calculation.

    returns array of strings
    '''
    mixedrep = loadtxt(findall('(?<=mixed representation\:).*', open(file).read()))[0::1]
    return mixedrep


def _twoCoscvelrep(file):
    '''
    Extract oscillator strengths in velocity representation from 2C escf calculation.

    returns array of strings
    '''
    velrep = loadtxt(findall('(?<=velocity representation\:).*', open(file).read()))[0::3]
    return velrep


def _twoCosclenrep(file):
    '''
    Extract oscillator strengths in length representation from 2C escf calculation.

    returns array of strings
    '''
    lenrep = loadtxt(findall('(?<=length representation\:).*', open(file).read()))[0::3]
    return lenrep


def _twoCoscmixrep(file):
    '''
    Extract oscillator strengths in mixed representation from 2C escf calculation.

    returns array of strings
    '''
    mixedrep = loadtxt(findall('(?<=mixed representation\:).*', open(file).read()))[0::2]
    return mixedrep


def _twoCtauvelrep(file):
    '''
    Extract radiative lifetimes in velocity representation from 2C escf calculation.

    returns array of strings
    '''
    velrep = loadtxt(findall('(?<=velocity representation\:).*', open(file).read()))[1::3]
    return velrep


def _twoCtaulenrep(file):
    '''
    Extract radiative lifetimes in length representation from 2C escf calculation.

    returns array of strings
    '''
    lenrep = loadtxt(findall('(?<=length representation\:).*', open(file).read()))[1::3]
    return lenrep


def _twoCtaumixrep(file):
    '''
    Extract radiative lifetimes in mixed representation from 2C escf calculation.

    returns array of strings
    '''
    mixedrep = loadtxt(findall('(?<=mixed representation\:).*', open(file).read()))[1::2]
    return mixedrep


def oscvelrep(file):
    if twoCcheck(file):
        velrep = _twoCoscvelrep(file)
    else:
        velrep = _oneCoscvelrep(file)
    return velrep


def osclenrep(file):
    if twoCcheck(file):
        lenrep = _twoCosclenrep(file)
    else:
        lenrep = _oneCosclenrep(file)
    return lenrep


def oscmixrep(file):
    if twoCcheck(file):
        mixrep = _twoCoscmixrep(file)
    else:
        mixrep = _oneCoscmixrep(file)
    return mixrep


def tauvelrep(file):
    if twoCcheck(file):
        velrep = _twoCtauvelrep(file)
    else:
        velrep = nan
    return velrep


def taulenrep(file):
    if twoCcheck(file):
        lenrep = _twoCtaulenrep(file)
    else:
        lenrep = nan
    return lenrep


def taumixrep(file):
    if twoCcheck(file):
        mixrep = _twoCtaumixrep(file)
    else:
        mixrep = nan
    return mixrep


def ricc2etenergies(file):
    raw = findall('\|      [0-9].[0-9][0-9][0-9][0-9][0-9] \|      [0-9].[0-9][0-9][0-9][0-9][0-9] \|', open(file).read())
    cleaned = [compile(r"\|").sub("", s) for s in raw]
    etenergies = eVtocm(loadtxt(cleaned)[:, 0])
    return etenergies


def ricc2etosclenrep(file):
    osc = loadtxt(findall('(?<=       oscillator strength \(length gauge\)   \:).*', open(file).read()))
    return osc


def ricc2etoscvelrep(file):
    osc = loadtxt(findall('(?<=       oscillator strength \(velocity gauge\) \:).*', open(file).read()))
    return osc


def ricc2etoscmixrep(file):
    osc = loadtxt(findall('(?<=       oscillator strength \(mixed gauge\)    \:).*', open(file).read()))
    return osc

