from numpy import loadtxt
from re import findall


def extractexcitations(file):
    excitations = loadtxt(findall('(?<=Excitation energy \/ cm\^\(-1\)\:).*', open(file).read()))
    return excitations


def extractvelrep(file):
    velrep = loadtxt(findall('(?<=Oscillator strength\:\n\n    velocity representation\:).*', open(file).read()))
    return velrep


def extractlenrep(file):
    lenrep = loadtxt(findall('(?<=length representation\:).*', open(file).read()))[0::3]
    return lenrep


def extractmixrep(file):
    mixedrep = loadtxt(findall('(?<=mixed representation\:).*', open(file).read()))[0::2]
    return mixedrep
