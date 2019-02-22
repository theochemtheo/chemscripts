import numpy as np
import re
from parseFUNCS import *


'''
parseORCA provides some additional functionality to, and is meant to work in concert with cclib
'''


def SOCME(file):
    """
    Parse an ORCA 4.1.0 file and recover SOC Matrix Elements from a CIS/TDA or RPA calculation

    returns a dictionary containing the following (complex valued) arrays with keys:
        x, y, z : x, y and z axis SOCME
        mag : normalized SOCME
        rows = T, columns = S
    """

    with open(file, 'r') as incoming:
        line = next(incoming)

        while line:
            reg_match = _ORCA4100_SOCME_reg(line)
            if reg_match.NROOTS_flag:
                nroots = int(line.split()[-1])
                line_count = (nroots + 1) * nroots
                SOCME_mag = np.zeros([nroots,nroots + 1], dtype=np.float_)
                SOCME_x = np.zeros([nroots,nroots + 1], dtype=np.complex_)
                SOCME_y = np.zeros([nroots,nroots + 1], dtype=np.complex_)
                SOCME_z = np.zeros([nroots,nroots + 1], dtype=np.complex_)
            if reg_match.SOCME_flag:
                # Skip some lines following match
                for N in range(5):
                    line = next(incoming, None)
                lines_done = 0
                while lines_done < line_count:
                    current_triplet = np.floor_divide(lines_done, nroots + 1)
                    current_singlet = np.mod(lines_done, nroots + 1)
                    soc_line = line.split()
                    socme_X = np.complex(np.float(soc_line[8]), np.float(soc_line[10]))
                    socme_Y = np.complex(np.float(soc_line[13]), np.float(soc_line[15]))
                    socme_Z = np.complex(np.float(soc_line[3]), np.float(soc_line[5]))
                    socme_mag = np.real(np.linalg.norm([socme_X, socme_Y, socme_Z]))
                    print('{}, {}: {}, {}, {}, {}'.format(current_triplet, current_singlet, socme_X, socme_Y, socme_Z, socme_mag))
                    SOCME_x[current_triplet, current_singlet] = socme_X
                    SOCME_y[current_triplet, current_singlet] = socme_Y
                    SOCME_z[current_triplet, current_singlet] = socme_Z
                    SOCME_mag[current_triplet, current_singlet] = socme_mag
                    line = next(incoming, None)
                    lines_done += 1

            # End of regex
            line = next(incoming, None)

    SOCME = {}
    SOCME['x'] = SOCME_x
    SOCME['y'] = SOCME_y
    SOCME['z'] = SOCME_z
    SOCME['mag'] = SOCME_mag

    return SOCME

class _ORCA4100_SOCME_reg:
    """ Regex for Gaussian 09 gradients """
    _reg_NROOTS_flag = re.compile(r'.*Number of roots to be determined.*')
    _reg_SOCME_flag = re.compile(r'.*CALCULATED SOCME BETWEEN TRIPLETS AND SINGLETS*')

    __slots__ = ['NROOTS_flag', 'SOCME_flag']

    def __init__(self, line):
        self.NROOTS_flag = self._reg_NROOTS_flag.match(line)
        self.SOCME_flag = self._reg_SOCME_flag.match(line)

