import numpy as np
import re
from parseFUNCS import *


'''
parseGAUSSIAN provides some additional functionality to, and is meant to work in concert with cclib
'''


def gradient(file):
    """
    Parse a Gaussian 09 file and recover gradient vector

    returns a numpy array containing the gradient vector
    """

    vector = np.array([]).reshape(0, 3)

    with open(file, 'r') as incoming:
        line = next(incoming)

        while line:
            reg_match = _G09_GRAD_reg(line)

            # Number of atoms needed later on
            if reg_match.SNO_flag:
                # skip 5 lines
                for N in range(5):
                    line = next(incoming, None)
                # Count lines until the end of the table is printed
                atomcount = 0
                while '------------' not in line:
                    atomcount += 1
                    line = next(incoming, None)

            # Figure out which coupling is being printed
            if reg_match.deriv_flag:
                # skip the flag line and next 2
                for N in range(3):
                    line = next(incoming, None)
                # Once at gradient itself, capture the next natom lines
                for i in range(atomcount):
                    vector = np.append(vector, [np.array(line.split()[2:], dtype=float)], axis=0)
                    line = next(incoming, None)

            # End of regex
            line = next(incoming, None)

    # Gradient = -Force
    GRADVEC = -vector

    return GRADVEC


class _G09_GRAD_reg:
    """ Regex for Gaussian 09 gradients """
    _reg_nuclei_flag = re.compile(r'.*Standard orientation:.*')
    _reg_deriv_flag = re.compile(r'.*Forces \(Hartrees/Bohr.*')

    __slots__ = ['SNO_flag', 'deriv_flag']

    def __init__(self, line):
        self.SNO_flag = self._reg_nuclei_flag.match(line)
        self.deriv_flag = self._reg_deriv_flag.match(line)
