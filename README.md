<!-- badges: start -->
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![GitHub license](https://img.shields.io/github/license/Naereen/StrapDown.js.svg)](/LICENSE)
[![CRAN version](https://www.r-pkg.org/badges/version/wmm)](https://cran.r-project.org/package=wmm)
[![CRAN downloads](https://cranlogs.r-pkg.org/badges/grand-total/wmm)](https://cran.r-project.org/package=wmm)
[![R-CMD-check](https://github.com/wfrierson/wmm/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/wfrierson/wmm/actions/workflows/R-CMD-check.yaml)
[![Coveralls test coverage](https://coveralls.io/repos/github/wfrierson/wmm/badge.svg)](https://coveralls.io/r/wfrierson/wmm?branch=main)
[![Status at rOpenSci Software Peer Review](https://badges.ropensci.org/522_status.svg)](https://github.com/ropensci/software-review/issues/522)
<!-- badges: end -->

# wmm
The [World Magnetic Model](https://www.ncei.noaa.gov/products/world-magnetic-model) (WMM)

The purpose of this package is to make accessible the magnetic field vector components from WMM. The supported date ranges for `wmm` run from 2000-01-01 to 2024-12-31. The magnetic field calculations across this time range agree with the official WMM test values to the precision provided by the authors (see the unit tests of this package for details). I will update this package for each new WMM version. For those that prefer a non-R solution, the authors of WMM provide free software to calculate magnetic field on the official website.

# Installation

Install from CRAN:

``` 
install.packages('wmm')
```

Or, install from GitHub:

``` r
devtools::install_github('wfrierson/wmm')

```

# Usage

In v1.1.2, the only exported function is `GetMagneticFieldWMM`, which returns the orthogonal vector components of the main magnetic field (in nT) and secular variation field (in nT/yr) as predicted by WMM. The magnetic field elements, _h_, _f_, _i_, and _d_ (as well as their secular variation) are returned as well.

Example usage: 

1. Calculate expected magnetic field components at a benchmark location, mid 2022. Using the default value for WMM version, this will use the WMM2020 coefficients.
``` r
wmm::GetMagneticFieldWMM(
  lon = 240,
  lat = -80,
  height = 1e5,
  time = 2022.5
)
# $x
# [1] 5814.966
#
# $y
# [1] 14802.97
#
# $z
# [1] -49755.31
#
# $xDot
# [1] 28.0382
#
# $yDot
# [1] 1.397062
#
# $zDot
# [1] 85.63095
#
# $h
# [1] 15904.14
#
# $f
# [1] 52235.36
#
# $i
# [1] -72.27367
#
# $d
# [1] 68.55389
#
# $hDot
# [1] 11.55182
#
# $fDot
# [1] -78.04815
#
# $iDot
# [1] 0.04066726
#
# $dDot
# [1] -0.09217566
```

2. Repeat the last calculation but apply it to 2017.5 and use the older coefficients from WMM2015 that were replaced by WMM2015v2. 

__Note: The WMM is intended to be predictive. By using an older set of coefficients, the returned values will reflect the older predictions. If users need a good model of the Earth's magnetic field prior to the current WMM, please see the latest [IGRF](https://www.ngdc.noaa.gov/IAGA/vmod/igrf.html) model, which is retroactively updated. The `wmmVersion` feature is intended for reproducibility purposes only.__
``` r
wmm::GetMagneticFieldWMM(
  lon = 240,
  lat = -80,
  height = 1e5,
  time = 2017.5,
  wmmVersion = 'WMM2015'
)
# $x
# [1] 5683.518
# 
# $y
# [1] 14808.85
# 
# $z
# [1] -50163.01
# 
# $xDot
# [1] 28.16496
# 
# $yDot
# [1] 6.941152
# 
# $zDot
# [1] 86.24356
# 
# $h
# [1] 15862.04
# 
# $f
# [1] 52611.14
# 
# $i
# [1] -72.45253
# 
# $d
# [1] 69.0036
# 
# $hDot
# [1] 16.57205
# 
# $fDot
# [1] -77.23403
# 
# $iDot
# [1] 0.04552524
# 
# $dDot
# [1] -0.08599694
```

# Citations

1. Chulliat, A., W. Brown, P. Alken, C. Beggan, M. Nair, G. Cox, A. Woods, 
S. Macmillan, B. Meyer and M. Paniccia, The US/UK World Magnetic Model for 
2020-2025: Technical Report, National Centers for Environmental Information, 
NOAA, doi: 10.25923/ytk1-yx35, 2020.

2. Chulliat, A., W. Brown, P. Alken, S. Macmillan, M. Nair, C. Beggan, A. Woods,
B. Hamilton, B. Meyer and R. Redmon, 2019, Out-of-Cycle Update of the
US/UK World Magnetic Model for 2015-2020: Technical Note, National
Centers for Environmental Information, NOAA. doi: 10.25921/xhr3-0t19

3. Chulliat, A., S. Macmillan, P. Alken, C. Beggan, M. Nair, B. Hamilton, A.
Woods, V. Ridley, S. Maus and A. Thomson, 2015, The US/UK World
Magnetic Model for 2015-2020: Technical Report, National Geophysical
Data Center, NOAA. doi: 10.7289/V5TB14V7

4. Maus, S., S. Macmillan, S. McLean, B. Hamilton, A. Thomson,
M. Nair, and C. Rollins, 2010, The US/UK World Magnetic Model
for 2010-2015, NOAA Technical Report NESDIS/NGDC.

5. McLean, S., S. Macmillan, S. Maus, V. Lesur, A.
Thomson, and D. Dater, December 2004, The
US/UK World Magnetic Model for 2005-2010,
NOAA Technical Report NESDIS/NGDC-1. 

6. Macmillian, S. and J. M. Quinn, 2000. 
“The Derivation of the World Magnetic Model 2000,” 
British Geological Survey Technical Report WM/00/17R.
