from numpy import pi, exp, log


def IRLorentzian(x, band, strength, FWHM):
    # Lorentzian lineshape function with IR prefactor
    # This is based on a document by Jens Spanget-Larsen, which may be found at dx.doi.org/10.13140/RG.2.1.4181.6160
    HWHM = FWHM / 2
    bandshape = ((100 / log(10) / pi * strength) / HWHM) / (1 + ((x - band) / HWHM)**2)
    return bandshape


def IRGaussian(x, band, strength, FWHM):
    # Gaussian lineshape function with IR prefactor
    HWHM = FWHM / 2
    bandshape = ((100 / log(10) / pi * strength) / HWHM) * exp(-2.772 * ((x - band) / FWHM)**2)
    return bandshape


def UVLorentzian(x, band, strength, FWHM):
    # Lorentzian lineshape function with UV prefactor
    # This is based on a document by Jens Spanget-Larsen, which may be found at dx.doi.org/10.13140/RG.2.1.4181.6160
    HWHM = FWHM / 2
    bandshape = (2.175e8) * (strength / FWHM) / (1 + ((x - band) / HWHM)**2)
    return bandshape


def UVGaussian(x, band, strength, FWHM):
    # Gaussian lineshape function with UV prefactor
    bandshape = (2.175e8) * (strength / FWHM) * exp(-2.772 * ((x - band) / FWHM)**2)
    return bandshape


def emGaussian(x, band, strength, FWHM):
    '''
    This provides a Gaussian lineshape for emission which is scaled in arbitrary units apart from
    the fact that the emission intensity should be cubically dependent on the transition energy
    This is due to the Einstein spontaneous emission coefficient, A.
    '''
    bandshape = x**3 * strength / FWHM * exp(-2.772 * ((x - band) / FWHM)**2)
    return bandshape
