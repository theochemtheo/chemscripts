import numpy as np
import re
from collections import OrderedDict
from parseFUNCS import *


'''

parseQCHEM provides additional functionality to, and is intended for use in concert with cclib

'''


def NACV(file):
    """
    Parse a Qchem 4.4 or 5.0 file and recover non-adiabatic coupling vectors

    returns two dictionaries, one containing the NACVs, h_ij (<Psi_i|d/dR|Psi_k>),
    NACVS['Mi-to-Mj'] = NACV between states i and j of multiplicity M
    these are off-diagonal elements of the force and thus have units hartrees bohr-1

    the other containing ancilliary data
    EXTRAS['calc_mult'] = multiplicity as S or T
    EXTRAS['calc_numstate'] = number of states between which NACVs calculated
    EXTRAS['calc_states'] = list of states between which couplings calculated
    EXTRAS['state_energies'] = list of state energies
    EXTRAS['deltaE_Mi-to-Mj'] = set delta E values in case of later need
    """
    NACVS = OrderedDict()
    EXTRAS = OrderedDict()
    EXTRAS['state_energies'] = []
    # defaults
    EXTRAS['calc_mult'] = 'S'
    with open(file, 'r') as incoming:
        line = next(incoming)

        while line:
            reg_match = _QC_NACV_reg(line)

            # If not a triplet, then singlet
            if reg_match.triplet:
                EXTRAS['calc_mult'] = 'T'

            # Number of states NACV is calculated for
            if reg_match.num_state:
                EXTRAS['calc_numstate'] = reg_match.num_state.group().split()[1]

            # The states between which NACVs are calculated
            if reg_match.states:
                line = next(incoming, None)
                line = next(incoming, None)
                state_string = re.search(r'^[^!]*', line)
                EXTRAS['calc_states'] = np.array([int(i) for i in state_string.group().split()])

            # Get the state energies so can divide derivative coupling with ETF by delta Eij later
            if reg_match.state_energy:
                energy = re.findall(r'[+-][0-9]*[.][0-9]+', line)[0]
                EXTRAS['state_energies'].append(float(energy))

            # Number of atoms needed later on
            if reg_match.SNO_flag:
                line = next(incoming, None)
                line = next(incoming, None)
                line = next(incoming, None)
                atomcount = 0
                while '------------' not in line:
                    atomcount += 1
                    line = next(incoming, None)
                EXTRAS['natom'] = atomcount

            # Figure out which coupling is being printed
            if reg_match.coupling_flag:
                state_i, state_j = re.findall('[0-9]+', line)
                NACV_key = '{0}{1}-to-{0}{2}'.format(EXTRAS['calc_mult'], state_i, state_j)
                vector = np.array([]).reshape(0, 3)
                # Loop through the until derivative coupling with ETF is printed
                while True:
                    line = next(incoming, None)
                    # Once at derivative coupling with ETF, skip through header then get next natoms lines and keep only the vector part
                    if re.search(r'.*CIS derivative coupling with ETF.*', line):
                        line = next(incoming, None)
                        line = next(incoming, None)
                        line = next(incoming, None)
                        for i in range(EXTRAS['natom']):
                            vector = np.append(vector, [np.array(line.split()[1:], dtype=float)], axis=0)
                            line = next(incoming, None)
                        break
                # Figure out delta Eij and keep
                delta_E = EXTRAS['state_energies'][int(state_j) - 1] - EXTRAS['state_energies'][int(state_i) - 1]
                EXTRAS['deltaE_{}'.format(NACV_key)] = delta_E
                # Multiply the DCwithETF by delta Eij to get NACV rather than derivative coupling
                vector = np.multiply(vector, delta_E)
                NACVS[NACV_key] = vector

            # End of regex
            line = next(incoming, None)

    return NACVS, EXTRAS


class _QC_NACV_reg:
    """ Regex for Qchem 4.4 & 5.0 NACVs """
    _reg_triplet = re.compile(r'cis_triplets.*true', re.IGNORECASE)
    _reg_numstate = re.compile(r'cis_der_numstate.*', re.IGNORECASE)
    _reg_states = re.compile(r'\$derivative_coupling', re.IGNORECASE)
    _reg_coupling_flag = re.compile(r'.*between states.*')
    _reg_state_energy = re.compile(r'.*Total energy for state.*')
    _reg_nuclei_flag = re.compile(r'.*Standard Nuclear Orientation.*')

    __slots__ = ['triplet', 'num_state', 'states', 'coupling_flag', 'state_energy', 'SNO_flag']

    def __init__(self, line):
        self.triplet = self._reg_triplet.match(line)
        self.num_state = self._reg_numstate.match(line)
        self.states = self._reg_states.match(line)
        self.coupling_flag = self._reg_coupling_flag.match(line)
        self.state_energy = self._reg_state_energy.match(line)
        self.SNO_flag = self._reg_nuclei_flag.match(line)


def gradient(file):
    """
    Parse a Qchem 4.4 or 5.0 file and recover gradient vector

    returns a numpy array containing the gradient vector
    """

    vector = np.array([]).reshape(0, 3)

    with open(file, 'r') as incoming:
        line = next(incoming)

        while line:
            reg_match = _QC_GRAD_reg(line)

            # Number of atoms needed later on
            if reg_match.SNO_flag:
                line = next(incoming, None)
                line = next(incoming, None)
                line = next(incoming, None)
                atomcount = 0
                while '------------' not in line:
                    atomcount += 1
                    line = next(incoming, None)

            # Figure out which coupling is being printed
            if reg_match.deriv_flag:
                # skip the flag line and next 3
                line = next(incoming, None)
                line = next(incoming, None)
                line = next(incoming, None)
                line = next(incoming, None)
                # Once at gradient itself, capture the next natom lines
                for i in range(atomcount):
                    vector = np.append(vector, [np.array(line.split()[1:], dtype=float)], axis=0)
                    line = next(incoming, None)

            # End of regex
            line = next(incoming, None)

    # For some reason, two gradients are printed. At this time not clear which is correct. Assuming the 2nd one.
    GRADVEC = vector[atomcount:]

    return GRADVEC


class _QC_GRAD_reg:
    """ Regex for Qchem 4.4 & 5.0 gradients """
    _reg_nuclei_flag = re.compile(r'.*Standard Nuclear Orientation.*')
    _reg_deriv_flag = re.compile(r'.*total gradient after adding PCM contribution.*')

    __slots__ = ['SNO_flag', 'deriv_flag']

    def __init__(self, line):
        self.SNO_flag = self._reg_nuclei_flag.match(line)
        self.deriv_flag = self._reg_deriv_flag.match(line)
