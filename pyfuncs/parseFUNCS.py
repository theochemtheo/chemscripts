import numpy as np


'''

parseFUNCS provides some generic functions used in the other parsing tools in this collection

'''

# globals
_bohrtoangstr = 0.529177210
_hartreetoeV = 27.211386
_hartreetocm = 2.194745313e5
_eVtocm = 8.065544e3
_daltontokg = 1.660539040e-27


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


def daltontokg(x):
    converted = x * _daltontokg
    return converted


def kgtodalton(x):
    converted = x / _daltontokg
    return converted


def atomwise_dot(vector1, vector2):
    """
    This function takes the dot product between two vectors row by row, e.g. for comparing two sets of molecular coordinates
    or for comparing a gradient with a vibrational mode

    returns a float
    """
    if len(vector1) != len(vector2):
        total_dot = np.nan
        raise ValueError('The two inputs are not the same length')
    else:
        total_dot = np.dot(np.ravel(vector1), np.ravel(vector2))
    return total_dot


def mass_weight(vector, atoms):
    """
    This function takes in a vector (i.e. NACV or gradient) and returns the mass-weighted version of it,
    where mass weighting does the following for each atom: mass weighted displacement = sqrt(mass) * displacement

    this reads in a vector and a list of atomic numbers
        the vector should be an array of cartesian displacements, one for each atom
        the list of atomic numbers should be a list of integers
    these must be in the same order, atom wise

    the mass-weighted vector is returned as a numpy array and has been normalised back to it's previous length
    """

    # Obtain masses from atom list
    masses = np.array([z_to_mass(Z) for Z in atoms])

    # Obtain the norm of the original vector
    orig_vec_norm = np.linalg.norm(vector, 'fro')

    # Generate the un-normalised mw_vector
    raw_mw_vector = np.array([]).reshape(0, 3)
    for row in range(len(vector)):
        thisatom = np.multiply(np.sqrt(masses[row]), vector[row, :])
        raw_mw_vector = np.append(raw_mw_vector, [thisatom], axis=0)

    # Normalise the mw vector to match the original norm
    raw_mw_vector_norm = np.linalg.norm(raw_mw_vector, 'fro')
    mw_vector = np.multiply(np.divide(raw_mw_vector, raw_mw_vector_norm), orig_vec_norm)

    return mw_vector


def unmass_weight(vector, atoms):
    """
    This function takes in a mass-weighted vector (i.e. NACV or gradient) and returns the unmass-weighted version of it,
    where mass weighting does the following for each atom: mass weighted displacement = sqrt(mass) * displacement

    this reads in a vector and a list of atomic numbers
        the vector should be an array of cartesian displacements, one for each atom
        the list of atomic numbers should be a list of integers
    these must be in the same order, atom wise

    the mass-weighted vector is returned as a numpy array and has been normalised back to it's previous length
    """

    # Obtain masses from atom list
    masses = np.array([z_to_mass(Z) for Z in atoms])

    # Obtain the norm of the original vector
    orig_mw_vec_norm = np.linalg.norm(vector, 'fro')

    # Generate the un-normalised mw_vector
    raw_umw_vector = np.array([]).reshape(0, 3)
    for row in range(len(vector)):
        thisatom = np.divide(vector[row, :], np.sqrt(masses[row]))
        raw_umw_vector = np.append(raw_umw_vector, [thisatom], axis=0)

    # Normalise the mw vector to match the original norm
    raw_umw_vector_norm = np.linalg.norm(raw_umw_vector, 'fro')
    umw_vector = np.multiply(np.divide(raw_umw_vector, raw_umw_vector_norm), orig_mw_vec_norm)

    return umw_vector


def z_to_mass(atomic_number):
    """
    This function takes in an atomic number and returns the standard atomic mass
    this was obtained from the CIAAW website on 2018/06/04
    www.ciaaw.or/atomic-weights.htm

    The abridged mass was used (accurate to 5 signifcant figures)
    For those atoms for which a range of masses is provided online, the midpoint was taken and rounded up.
    For atoms with no stable isotopes, the mass of the most stable isotope was used as obtained from webelements
    www.webelements.com

    returns a float
    """

    if atomic_number == 1:
        atomic_mass = 1.0080
    elif atomic_number == 2:
        atomic_mass = 4.0026
    elif atomic_number == 3:
        atomic_mass = 6.9675
    elif atomic_number == 4:
        atomic_mass = 9.0122
    elif atomic_number == 5:
        atomic_mass = 10.8135
    elif atomic_number == 6:
        atomic_mass = 12.0105
    elif atomic_number == 7:
        atomic_mass = 14.007
    elif atomic_number == 8:
        atomic_mass = 16.000
    elif atomic_number == 9:
        atomic_mass = 18.998
    elif atomic_number == 10:
        atomic_mass = 20.180
    elif atomic_number == 11:
        atomic_mass = 22.990
    elif atomic_number == 12:
        atomic_mass = 24.306
    elif atomic_number == 13:
        atomic_mass = 26.982
    elif atomic_number == 14:
        atomic_mass = 28.085
    elif atomic_number == 15:
        atomic_mass = 30.974
    elif atomic_number == 16:
        atomic_mass = 32.068
    elif atomic_number == 17:
        atomic_mass = 35.452
    elif atomic_number == 18:
        atomic_mass = 39.948
    elif atomic_number == 19:
        atomic_mass = 39.098
    elif atomic_number == 20:
        atomic_mass = 40.078
    elif atomic_number == 21:
        atomic_mass = 44.956
    elif atomic_number == 22:
        atomic_mass = 47.867
    elif atomic_number == 23:
        atomic_mass = 50.942
    elif atomic_number == 24:
        atomic_mass = 51.996
    elif atomic_number == 25:
        atomic_mass = 54.938
    elif atomic_number == 26:
        atomic_mass = 55.845
    elif atomic_number == 27:
        atomic_mass = 58.933
    elif atomic_number == 28:
        atomic_mass = 58.693
    elif atomic_number == 29:
        atomic_mass = 63.546
    elif atomic_number == 30:
        atomic_mass = 65.380
    elif atomic_number == 31:
        atomic_mass = 69.723
    elif atomic_number == 32:
        atomic_mass = 72.630
    elif atomic_number == 33:
        atomic_mass = 74.922
    elif atomic_number == 34:
        atomic_mass = 78.971
    elif atomic_number == 35:
        atomic_mass = 79.904
    elif atomic_number == 36:
        atomic_mass = 83.798
    elif atomic_number == 37:
        atomic_mass = 85.468
    elif atomic_number == 38:
        atomic_mass = 87.620
    elif atomic_number == 39:
        atomic_mass = 88.906
    elif atomic_number == 40:
        atomic_mass = 91.224
    elif atomic_number == 41:
        atomic_mass = 92.906
    elif atomic_number == 42:
        atomic_mass = 95.950
    elif atomic_number == 43:
        atomic_mass = 98.000
    elif atomic_number == 44:
        atomic_mass = 101.07
    elif atomic_number == 45:
        atomic_mass = 102.91
    elif atomic_number == 46:
        atomic_mass = 106.42
    elif atomic_number == 47:
        atomic_mass = 107.87
    elif atomic_number == 48:
        atomic_mass = 112.41
    elif atomic_number == 49:
        atomic_mass = 114.82
    elif atomic_number == 50:
        atomic_mass = 118.71
    elif atomic_number == 51:
        atomic_mass = 121.76
    elif atomic_number == 52:
        atomic_mass = 127.60
    elif atomic_number == 53:
        atomic_mass = 126.90
    elif atomic_number == 54:
        atomic_mass = 131.29
    elif atomic_number == 55:
        atomic_mass = 132.91
    elif atomic_number == 56:
        atomic_mass = 137.33
    elif atomic_number == 57:
        atomic_mass = 138.91
    elif atomic_number == 58:
        atomic_mass = 140.12
    elif atomic_number == 59:
        atomic_mass = 140.91
    elif atomic_number == 60:
        atomic_mass = 144.24
    elif atomic_number == 61:
        atomic_mass = 144.91
    elif atomic_number == 62:
        atomic_mass = 150.36
    elif atomic_number == 63:
        atomic_mass = 151.96
    elif atomic_number == 64:
        atomic_mass = 157.25
    elif atomic_number == 65:
        atomic_mass = 158.93
    elif atomic_number == 66:
        atomic_mass = 162.50
    elif atomic_number == 67:
        atomic_mass = 164.93
    elif atomic_number == 68:
        atomic_mass = 167.26
    elif atomic_number == 69:
        atomic_mass = 168.93
    elif atomic_number == 70:
        atomic_mass = 173.05
    elif atomic_number == 71:
        atomic_mass = 174.97
    elif atomic_number == 72:
        atomic_mass = 178.49
    elif atomic_number == 73:
        atomic_mass = 180.95
    elif atomic_number == 74:
        atomic_mass = 183.84
    elif atomic_number == 75:
        atomic_mass = 186.21
    elif atomic_number == 76:
        atomic_mass = 190.23
    elif atomic_number == 77:
        atomic_mass = 192.22
    elif atomic_number == 78:
        atomic_mass = 195.08
    elif atomic_number == 79:
        atomic_mass = 196.97
    elif atomic_number == 80:
        atomic_mass = 200.59
    elif atomic_number == 81:
        atomic_mass = 204.39
    elif atomic_number == 82:
        atomic_mass = 207.20
    elif atomic_number == 83:
        atomic_mass = 208.98
    elif atomic_number == 84:
        atomic_mass = 208.98
    elif atomic_number == 85:
        atomic_mass = 209.99
    elif atomic_number == 86:
        atomic_mass = 222.02
    elif atomic_number == 87:
        atomic_mass = 223.02
    elif atomic_number == 88:
        atomic_mass = 226.03
    elif atomic_number == 89:
        atomic_mass = 227.03
    elif atomic_number == 90:
        atomic_mass = 232.04
    elif atomic_number == 91:
        atomic_mass = 231.04
    elif atomic_number == 92:
        atomic_mass = 238.03
    elif atomic_number == 93:
        atomic_mass = 237.05
    elif atomic_number == 94:
        atomic_mass = 244.06
    elif atomic_number == 95:
        atomic_mass = 243.06
    elif atomic_number == 96:
        atomic_mass = 247.07
    elif atomic_number == 97:
        atomic_mass = 247.07
    elif atomic_number == 98:
        atomic_mass = 251.08
    elif atomic_number == 99:
        atomic_mass = 252.08
    elif atomic_number == 100:
        atomic_mass = 257.10
    elif atomic_number == 101:
        atomic_mass = 258.10
    elif atomic_number == 102:
        atomic_mass = 259.10
    elif atomic_number == 103:
        atomic_mass = 262.11
    elif atomic_number == 104:
        atomic_mass = 267.12
    elif atomic_number == 105:
        atomic_mass = 270.13
    elif atomic_number == 106:
        atomic_mass = 269.13
    elif atomic_number == 107:
        atomic_mass = 270.13
    elif atomic_number == 108:
        atomic_mass = 269.13
    elif atomic_number == 109:
        atomic_mass = 278.16
    elif atomic_number == 110:
        atomic_mass = 281.17
    elif atomic_number == 111:
        atomic_mass = 281.17
    elif atomic_number == 112:
        atomic_mass = 285.18
    elif atomic_number == 113:
        atomic_mass = 286.18
    elif atomic_number == 114:
        atomic_mass = 289.19
    elif atomic_number == 115:
        atomic_mass = 289.20
    elif atomic_number == 116:
        atomic_mass = 293.20
    elif atomic_number == 117:
        atomic_mass = 293.21
    elif atomic_number == 118:
        atomic_mass = 294.21
    else:
        print('Error reading atomic number: {}'.format(atomic_number))
        exit()
    return atomic_mass


def osc_to_tdm(f, Ei, Ef):
    '''
    Convert oscillator strength to transition dipole moment

    requires initial and final energies in hartrees

    returns tdm in atomic units (e * bohr) as float
    '''
    if f == 0:
        tdm = 0
    else:
        tdm = np.sqrt((3 * f) / (2 * (Ef - Ei)))
    return tdm
