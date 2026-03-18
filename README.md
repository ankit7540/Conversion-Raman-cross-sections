# Conversion-Raman-cross-sections
Set of python functions to aid in converting Raman cross-section values including unit conversion.
Probably the most useful function is to get the Raman cross-section at different wavelengths, useful for comparison of literature data.


## Available functions

```
      RCS_interpolate_to_wavelength( RCS, freq , lambda_org, lambda_interp)

      localField_corr(n_exc, n_sc)

      differential_to_total_RCS(dep_ratio, diffRCS)

      convert_RCS_au_to_cm6(RCS)

      convert_RCS_au_to_SI(RCS)

      convert_RCS_au_to_Angstrom(RCS)
```


## References

Specific references are defined in the body of the functions.
