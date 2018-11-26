from numpy import pi, exp, log, sqrt
import parseFUNCS as pF


_daltontokg = 1.660539040e-27
_kBoltzmann = 1.38064852e-23
_hPlanck = 6.62607004e-34
_cSpeedOfLight = 299792458


def _IRprefac(oscillator_strength, FWHM):
    HWHM = FWHM / 2
    prefac = ((100 / log(10) / pi * oscillator_strength) / HWHM)
    return prefac


def _Ramanprefac(transition_energy, activity, laser_excitation, temperature):
    '''
    Prefactor for non-resonant Raman scattering intensity
    This is based on eqn. 1 of V. Krishnakumar et al., J. Mol. Struct. 2004, 702, 9
    dx.doi.org/10.1016/j.molstruc.2004.06.004
    '''
    prefac_num = (100 * (1e7 / laser_excitation - transition_energy))**4 * (1e-40 * activity / _daltontokg)
    prefac_denom = (100 * transition_energy) * (1 - exp(-(_hPlanck * _cSpeedOfLight * 100 * transition_energy) / (_kBoltzmann * temperature)))
    prefac = prefac_num / prefac_denom
    return prefac


def _STA_resonance_raman_prefac(transition_energy, tdm, gradient):
    '''
    Equation 8 from Kane & Jansen, J. Phys. Chem. C, 2010, 114, 5541 dx.doi.org/10.1021/jp906152q
    "Short time approximation" or ES-gradient approximation

    dimensionless normal mode coordinate q = sqrt(h-bar / frequency) * Normalmode
    displacement from equation 6 of same paper (IMDHO model)
    '''
    normal_mode_to_dimless_factor = sqrt(1 / pF.cmtohartree(transition_energy))
    dimless_gradient = gradient / normal_mode_to_dimless_factor
    displacement = dimless_gradient / -pF.cmtohartree(transition_energy)
    scattering_factor = 12 * tdm**4 * transition_energy**2 * displacement**2
    return scattering_factor


def _UVprefac(oscillator_strength, FWHM):
    '''
    This is based on a document by Jens Spanget-Larsen, which may be found at dx.doi.org/10.13140/RG.2.1.4181.6160
    '''
    prefac = (2.175e8) * (oscillator_strength / FWHM)
    return prefac


def _PhotoAbsCrossSec(x, transition_energy, oscillator_strength):
    '''
    This function implements the photoabsorption cross section, from which other observables can be derived.
    The incoming values must be in hartree atomic units

    This is based on equation 87 in Crespo-Otero & Barbatti's Chem. Rev., 2018, 118, 7026 dx.doi/10.1021/acs.chemrev.7b00577
    '''
    prefac = (2 * pi**2) / (137 * x)
    cross_section = prefac * transition_energy * oscillator_strength
    return cross_section


def _emprefac(x, oscillator_strength, FWHM):
    prefac = x**3 * oscillator_strength / FWHM
    return prefac


def _Lorentzian(x, transition_energy, FWHM):
    HWHM = FWHM / 2
    probability = 1 / (1 + ((x - transition_energy) / HWHM)**2)
    return probability


def _Gaussian(x, transition_energy, FWHM):
    probability = exp(-2.772 * ((x - transition_energy) / FWHM)**2)
    return probability


def IRLorentzian(x, transition_energy, oscillator_strength, FWHM):
    '''
    Lorentzian lineshape function with IR prefactor
    This is based on a document by Jens Spanget-Larsen, which may be found at dx.doi.org/10.13140/RG.2.1.4181.6160
    '''
    epsilon = _IRprefac(oscillator_strength, FWHM) * _Lorentzian(x, transition_energy, FWHM)
    return epsilon


def IRGaussian(x, transition_energy, oscillator_strength, FWHM):
    '''
    Gaussian lineshape function with IR prefactor
    '''
    epsilon = _IRprefac(oscillator_strength, FWHM) * _Gaussian(x, transition_energy, FWHM)
    return epsilon


def RamanLorentzian(x, transition_energy, activity, laser_excitation, temperature, FWHM):
    '''
    Lorentzian lineshape function with Raman prefactor
    '''
    intensity = _Ramanprefac(transition_energy, activity, laser_excitation, temperature) * _Lorentzian(x, transition_energy, FWHM)
    return intensity


def RamanGaussian(x, transition_energy, activity, laser_excitation, temperature, FWHM):
    '''
    Gaussian lineshape function with Raman prefactor
    '''
    intensity = _Ramanprefac(transition_energy, activity, laser_excitation, temperature) * _Gaussian(x, transition_energy, FWHM)
    return intensity


def UVLorentzian(x, transition_energy, oscillator_strength, FWHM):
    '''
    Lorentzian lineshape function with UV prefactor
    '''
    epsilon = _UVprefac(oscillator_strength, FWHM) * _Lorentzian(x, transition_energy, FWHM)
    return epsilon


def UVGaussian(x, transition_energy, oscillator_strength, FWHM):
    '''
    Gaussian lineshape function with UV prefactor
    '''
    epsilon = _UVprefac(oscillator_strength, FWHM) * _Gaussian(x, transition_energy, FWHM)
    return epsilon


def emLorentzian(x, transition_energy, oscillator_strength, FWHM):
    '''
    This provides a Lorentzian lineshape for emission which is scaled in arbitrary units apart from
    the fact that the emission intensity should be cubically dependent on the transition energy
    This is due to the Einstein spontaneous emission coefficient, A.
    '''
    epsilon = _emprefac(x, oscillator_strength, FWHM) * _Lorentzian(x, transition_energy, FWHM)
    return epsilon


def emGaussian(x, transition_energy, oscillator_strength, FWHM):
    '''
    This provides a Gaussian lineshape for emission which is scaled in arbitrary units apart from
    the fact that the emission intensity should be cubically dependent on the transition energy
    This is due to the Einstein spontaneous emission coefficient, A.
    '''
    epsilon = _emprefac(x, oscillator_strength, FWHM) * _Gaussian(x, transition_energy, FWHM)
    return epsilon


def PhotoAbsCrossSectionGaussian(x, transition_energy, oscillator_strength, FWHM):
    '''
    This provides a Gaussian lineshape for photoabsorption cross section. Note that the incoming values must be in atomic units

    returns cross section (float)
    '''
    cross_section = _PhotoAbsCrossSec(x, transition_energy, oscillator_strength) * _Gaussian(x, transition_energy, FWHM)
    return cross_section


def PhotoAbsCrossSectionLorentzian(x, transition_energy, oscillator_strength, FWHM):
    '''
    This provides a Lorentzian lineshape for photoabsorption cross section. Note that the incoming values must be in atomic units

    returns cross section (float)
    '''
    cross_section = _PhotoAbsCrossSec(x, transition_energy, oscillator_strength) * _Lorentzian(x, transition_energy, FWHM)
    return cross_section


def STA_RR_Gaussian(x, transition_energy, tdm, gradient, FWHM):
    '''
    Resonance Raman within the short time approximation a.k.a. ES-gradient approximation broadened with a Gaussian

    returns intensity (float)
    '''
    intensity = _STA_resonance_raman_prefac(transition_energy, tdm) * _Gaussian(x, transition_energy, FWHM)
    return intensity


def STA_RR_Lorentzian(x, transition_energy, tdm, gradient, FWHM):
    '''
    Resonance Raman within the short time approximation a.k.a. ES-gradient approximation broadened with a Lorentzian

    returns intensity (float)
    '''
    intensity = _STA_resonance_raman_prefac(transition_energy, tdm) * _Lorentzian(x, transition_energy, FWHM)
    return intensity
