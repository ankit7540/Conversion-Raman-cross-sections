# Conversion-Raman-cross-sections

Set of python functions to aid in converting Raman cross-section values measured at certain wavelength to an approximate value at another wavelength. This conversion accounts only for the frequency factor and the effect of pre-resonance amplification is not included.

Additionally, functions for unit conversion of RCS values are included.

See function definitions for usage.


## Available functions

```
      RCS_interpolate_to_wavelength( RCS, freq , lambda_org, lambda_interp)

      RCS_interpolate_to_wavelength_PhotonCounting( RCS, freq , lambda_org, lambda_interp)

      RCS_interpolate_to_wavelength_power( RCS, freq , lambda_org, lambda_interp)

      localField_corr(n_exc, n_sc)

      differential_to_total_RCS(dep_ratio, diffRCS)

      convert_RCS_au_to_cm6(RCS)

      convert_RCS_au_to_SI(RCS)

      convert_RCS_au_to_Angstrom(RCS)
```

## Requirements

`Python3`


## References

Specific references are defined in the body of the functions.


## Used in 

 - [Accurate determination of absolute Raman cross-sections of carbon dioxide using H2 as reference](https://doi.org/10.1021/acs.jpca.6c00827)
