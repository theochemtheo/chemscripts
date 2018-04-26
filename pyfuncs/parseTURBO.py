from numpy import loadtxt, nan, array, append, float
from re import findall, compile


# globals
_bohrtoangstr = 0.529177210
_hartreetoeV = 27.211386
_hartreetocm = 2.194745313e5
_eVtocm = 8.065544e3


def atomcoords(file):
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
    gs = loadtxt(findall('(?<=Total energy\:).*', open(file).read()))[0]
    return gs


def esenergies(file):
    es = loadtxt(findall('(?<=Total energy\:).*', open(file).read()))[1:]
    return es


def excitations(file):
    excitations = loadtxt(findall('(?<=Excitation energy \/ cm\^\(-1\)\:).*', open(file).read()))
    return excitations


def twoCcheck(file):
    istwoC = False
    with open(file, 'r') as f:
        for line in f.readlines():
            if 'Two-component modus switched on ! ' in line:
                istwoC = True
    return istwoC


def _oneCoscvelrep(file):
    velrep = loadtxt(findall('(?<=velocity representation\:).*', open(file).read()))[0::2]
    return velrep


def _oneCosclenrep(file):
    lenrep = loadtxt(findall('(?<=length representation\:).*', open(file).read()))[0::2]
    return lenrep


def _oneCoscmixrep(file):
    mixedrep = loadtxt(findall('(?<=mixed representation\:).*', open(file).read()))[0::1]
    return mixedrep


def _twoCoscvelrep(file):
    velrep = loadtxt(findall('(?<=velocity representation\:).*', open(file).read()))[0::3]
    return velrep


def _twoCosclenrep(file):
    lenrep = loadtxt(findall('(?<=length representation\:).*', open(file).read()))[0::3]
    return lenrep


def _twoCoscmixrep(file):
    mixedrep = loadtxt(findall('(?<=mixed representation\:).*', open(file).read()))[0::2]
    return mixedrep


def _twoCtaulenrep(file):
    lenrep = loadtxt(findall('(?<=length representation\:).*', open(file).read()))[1::3]
    return lenrep


def _twoCtaumixrep(file):
    mixedrep = loadtxt(findall('(?<=mixed representation\:).*', open(file).read()))[1::2]
    return mixedrep


def _twoCtauvelrep(file):
    velrep = loadtxt(findall('(?<=velocity representation\:).*', open(file).read()))[1::3]
    return velrep


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


def cmtoeV(x):
    converted = x / _eVtocm
    return converted


def eVtocm(x):
    converted = x * _eVtocm
    return converted


def hartreetocm(x):
    converted = x * _hartreetocm
    return converted


def cmtohartree(x):
    converted = x / _hartreetocm
    return converted


def hartreetoeV(x):
    converted = x * _hartreetoeV
    return converted


def eVtohartree(x):
    converted = x / _hartreetoeV
    return converted


def bohrtoangstr(x):
    converted = x * _bohrtoangstr
    return converted


def angstrtobohr(x):
    converted = x / _bohrtoangstr
    return converted
