from numpy import pi, exp, log


_daltontokg = 1.660539040e-27
_kBoltzmann = 1.38064852e-23
_hPlanck = 6.62607004e-34
_cSpeedOfLight = 299792458


def _IRprefac(strength, FWHM):
    HWHM = FWHM / 2
    prefac = ((100 / log(10) / pi * strength) / HWHM)
    return prefac


def _Ramanprefac(band, activity, excitation, temperature):
    '''
    Prefactor for non-resonant Raman scattering intensity
    This is based on eqn. 1 of V. Krishnakumar et al., J. Mol. Struct. 2004, 702, 9
    dx.doi.org/10.1016/j.molstruc.2004.06.004
    '''
    prefac_num = (100 * (1e7 / excitation - band))**4 * (1e-40 * activity / _daltontokg)
    prefac_denom = (100 * band) * (1 - exp(-(_hPlanck * _cSpeedOfLight * 100 * band) / (_kBoltzmann * temperature)))
    prefac = prefac_num / prefac_denom
    return prefac


def _UVprefac(strength, FWHM):
    '''
    This is based on a document by Jens Spanget-Larsen, which may be found at dx.doi.org/10.13140/RG.2.1.4181.6160
    '''
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
    '''
    Lorentzian lineshape function with IR prefactor
    This is based on a document by Jens Spanget-Larsen, which may be found at dx.doi.org/10.13140/RG.2.1.4181.6160
    '''
    epsilon = _IRprefac(strength, FWHM) * _Lorentzian(x, band, FWHM)
    return epsilon


def IRGaussian(x, band, strength, FWHM):
    '''
    Gaussian lineshape function with IR prefactor
    '''
    epsilon = _IRprefac(strength, FWHM) * _Gaussian(x, band, FWHM)
    return epsilon


def RamanLorentzian(x, band, activity, excitation, temperature, FWHM):
    '''
    Lorentzian lineshape function with Raman prefactor
    '''
    intensity = _Ramanprefac(band, activity, excitation, temperature) * _Lorentzian(x, band, FWHM)
    return intensity


def RamanGaussian(x, band, activity, excitation, temperature, FWHM):
    '''
    Gaussian lineshape function with Raman prefactor
    '''
    intensity = _Ramanprefac(band, activity, excitation, temperature) * _Gaussian(x, band, FWHM)
    return intensity


def UVLorentzian(x, band, strength, FWHM):
    '''
    Lorentzian lineshape function with UV prefactor
    '''
    epsilon = _UVprefac(strength, FWHM) * _Lorentzian(x, band, FWHM)
    return epsilon


def UVGaussian(x, band, strength, FWHM):
    '''
    Gaussian lineshape function with UV prefactor
    '''
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
