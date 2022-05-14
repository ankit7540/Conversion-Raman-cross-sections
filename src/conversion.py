


###################################################################

def RCS_interpolate_to_wavelength( RCS, freq , lambda_org, lambda_interp):
    '''
    To convert a given Raman cross-section measured at a wavelength value to
    approximate value for a different excitation laser

    Input params :
    RCS = value of Raman cross-section , for example, 1.0e-30 or
                                         scaled value like 1.0
    freq = vibrational frequency of the specific transition (in relative wavenumbers)
            for example, 992.3
    lambda_org = wavelength in nm, for the excitation laser
    lambda_interp = wavelength in nm, for the laser at which we want to approximate the Raman Cross-section
    '''
    # --------------
    # compute freq_factor for original wavelength
    wavenum_laser = 1e7/lambda_org
    wavenum_sc = wavenum_laser - freq

    #print (wavenum_laser, wavenum_sc)
    freq_factor = wavenum_laser * (wavenum_sc**3)

    # --------------
    # compute freq_factor for the new wavelength
    wavenum_interp = 1e7/lambda_interp
    wavenum_interp_sc = wavenum_interp - freq

    freq_factor_interp = wavenum_interp * (wavenum_interp_sc**3)

    print (" Scaling factor : ", ((1 / freq_factor) * freq_factor_interp) )

    return RCS * (1/freq_factor) * freq_factor_interp

###################################################################

def localField_corr(n_exc, n_sc):
    '''To compute the local-field correction factor for ref. indices at given specific
    excitation laser and the scattering wavelength

    Input Params :
    n_exc = ref. index of the material at excitation laser. for example, 1.33
    n_sc  = ref. index of the material at the scattering wavelength, for example, 1.4

    Ref.
        Nestor & Lippincott :  https://doi.org/10.1002/jrs.1250010309
        Eckhardt & Wagner   :  https://doi.org/10.1016/0022-2852(66)90262-1
        Abe, Wakayama & Ito :   https://doi.org/10.1002/jrs.1250060109

    '''

    val_exc = ((n_exc**2)+2)
    val_sc = ((n_sc**2)+2)
    return (val_sc**2)*(val_exc**2)*(n_sc/n_exc)*(1/81)


###################################################################


def differential_to_total_RCS(dep_ratio, diffRCS):
    '''
    To convert differential Raman cross-section to corresponding total Raman
    cross-section value using the depolarization ratio

    Input Params:
    dep_ratio  = depolarization ratio
    diffRCS    = differential Raman cross-section value

    Ref.
            Trulson and Mathies  :  https://doi.org/10.1063/1.450415
            W. R. Hess, H. Hacker, H. W. Schrotter, and J. Brandmiiller, Z. Angew.
              Phys. 27, 233 (1969)


    '''

    pi = 3.14159265358979
    return (8* pi/3) * ((1+2*dep_ratio)/(1+dep_ratio)) * (diffRCS)

###################################################################
