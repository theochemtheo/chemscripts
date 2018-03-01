from numpy import loadtxt
from re import findall


def excitations(file):
    excitations = loadtxt(findall('(?<=Excitation energy \/ cm\^\(-1\)\:).*', open(file).read()))
    return excitations


def oscvelrep(file):
    velrep = loadtxt(findall('(?<=velocity representation\:).*', open(file).read()))[0::3]
    return velrep


def osclenrep(file):
    lenrep = loadtxt(findall('(?<=length representation\:).*', open(file).read()))[0::3]
    return lenrep


def oscmixrep(file):
    mixedrep = loadtxt(findall('(?<=mixed representation\:).*', open(file).read()))[0::2]
    return mixedrep


def taulenrep(file):
    lenrep = loadtxt(findall('(?<=length representation\:).*', open(file).read()))[1::3]
    return lenrep


def taumixrep(file):
    mixedrep = loadtxt(findall('(?<=mixed representation\:).*', open(file).read()))[1::2]
    return mixedrep


def tauvelrep(file):
    velrep = loadtxt(findall('(?<=velocity representation\:).*', open(file).read()))[1::3]
    return velrep
