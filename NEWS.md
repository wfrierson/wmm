# wmm 1.1.3
1. Added `inst/extdata/coefficients/` to store legacy coefficient files, though these are included in `.gitignore`.
2. Added WMM2025 coefficient path to `data-raw/creatingSysData.R`.
3. Added WMM2025 test values to `inst/extdata/WMMTestValues.csv`.
4. Recalculated `R/sysdata.rda` to incorporate WMM2025 coefficients.
5. Removed legacy roxygen syntax in `R/wmm.R`.
6. In `R/functions_wmm.R`, updated roxygen labels and `wmmVersion` @param descriptions to include 'WMM2025', where appropriate.
7. In `R/functions_wmm.R`, included 'WMM2000' in `wmmVersion` @param description of `.CalculateMagneticField`.
8. In `R/functions_wmm.R`, updated @examples for `GetMagneticFieldWMM` to use an example from the official WMM2025 test values.
9. In `R/functions_misc.R`, added date ranges for WMM2025 and WMM2000 coefficients in `.DeriveVersionInfo`.
10. In `R/functions_misc.R`, included 'WMM2025' in `wmmVersion` @param description of `.CheckVersionWMM`.
11. In `R/functions_misc.R`, updated URL for 'Performance Specifications WMM' in `.CheckBlackoutZone`.
12. In `R/functions_coefficients.R`, included 'WMM2025' in `wmmVersion` @param description of `.CalculateGaussCoef`.
13. Updated documentation in manual.
14. Updated `tests/testthat.R` to work in vscode.
15. Updated `vignettes/` to reference WMM2025, where appropriate.
16. In `DESCRIPTION`:
i. updated supported date ranges
ii. included WMM2025 citation
iii. updated versions within `Depends` and `Suggests`
17. In `README.md`:
i. updated supported date ranges
ii. included WMM2025 citation
18. Updated `renv.lock` to use latest package versions and R 4.3.3.

# wmm 1.1.2
1. Added HTML vignette.
2. Added @noRd flag for internal functions.
3. Added docs/CONTRIBUTING.md
4. Added codemeta.json file.

# wmm 1.1.1
1. Removed WMM coefficient files and URLs.
2. Adding note to README re: WMM output as predictions and referenced IGRF.

# wmm 1.1.0
1. Updated Gauss coefficients for WMM2020.
2. Significantly improved numerical stability by using closed-form equation for associated Legendre polynomials, instead of recursion.
3. Improved speed by ~7x by using multidimensional arrays.
4. Per WMM2020, this version now displays a blackout zone warning when the horizontal intensity is below 6000 nT.
5. In addition to orthogonal magnetic field components, `GetMagneticFieldWMM` returns the magnetic field elements _h_, _f_, _i_, & _d_ as well as their secular variation.

# wmm 1.0.0
First release of wmm package. Contains 1 exported function, `GetMagneticFieldWMM`, which calculates magnetic field.
