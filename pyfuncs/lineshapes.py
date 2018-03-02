from numpy import pi, exp, log


def _IRprefac(strength, FWHM):
    HWHM = FWHM / 2
    prefac = ((100 / log(10) / pi * strength) / HWHM)
    return prefac


def _UVprefac(strength, FWHM):
    prefac = (2.175e8) * (strength / FWHM)
    return prefac


def _emprefac(x, strength, FWHM):
    prefac = x**3 * strength / FWHM
    return prefac


def _Lorentzian(x, band, FWHM):
    HWHM = FWHM / 2
    probability = 1 / (1 + ((x - band) / HWHM)**2)
    return probability


def _Gaussian(x, band, FWHM):
    probability = exp(-2.772 * ((x - band) / FWHM)**2)
    return probability


def IRLorentzian(x, band, strength, FWHM):
    # Lorentzian lineshape function with IR prefactor
    # This is based on a document by Jens Spanget-Larsen, which may be found at dx.doi.org/10.13140/RG.2.1.4181.6160
    epsilon = _IRprefac(strength, FWHM) * _Lorentzian(x, band, FWHM)
    return epsilon


def IRGaussian(x, band, strength, FWHM):
    # Gaussian lineshape function with IR prefactor
    epsilon = _IRprefac(strength, FWHM) * _Gaussian(x, band, FWHM)
    return epsilon


def UVLorentzian(x, band, strength, FWHM):
    # Lorentzian lineshape function with UV prefactor
    # This is based on a document by Jens Spanget-Larsen, which may be found at dx.doi.org/10.13140/RG.2.1.4181.6160
    epsilon = _UVprefac(strength, FWHM) * _Lorentzian(x, band, FWHM)
    return epsilon


def UVGaussian(x, band, strength, FWHM):
    # Gaussian lineshape function with UV prefactor
    epsilon = _UVprefac(strength, FWHM) * _Gaussian(x, band, FWHM)
    return epsilon


def emLorentzian(x, band, strength, FWHM):
    '''
    This provides a Lorentzian lineshape for emission which is scaled in arbitrary units apart from
    the fact that the emission intensity should be cubically dependent on the transition energy
    This is due to the Einstein spontaneous emission coefficient, A.
    '''
    epsilon = _emprefac(x, strength, FWHM) * _Lorentzian(x, band, FWHM)
    return epsilon


def emGaussian(x, band, strength, FWHM):
    '''
    This provides a Gaussian lineshape for emission which is scaled in arbitrary units apart from
    the fact that the emission intensity should be cubically dependent on the transition energy
    This is due to the Einstein spontaneous emission coefficient, A.
    '''
    epsilon = _emprefac(x, strength, FWHM) * _Gaussian(x, band, FWHM)
    return epsilon
