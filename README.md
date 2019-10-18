<!-- badges: start -->
[![Travis build status](https://travis-ci.org/wfrierson/wmm.svg?branch=master)](https://travis-ci.org/wfrierson/wmm)
[![GitHub license](https://img.shields.io/github/license/Naereen/StrapDown.js.svg)](/LICENSE)
[![Coverage Status](https://coveralls.io/repos/github/wfrierson/wmm/badge.svg?branch=master)](https://coveralls.io/github/wfrierson/wmm?branch=master)
<!-- badges: end -->

# wmm
Unofficial R implementation of the [World Magnetic Model](https://www.ngdc.noaa.gov/geomag/WMM/DoDWMM.shtml) (WMM). This is a work in progress.

The purpose of this project is to program WMM entirely in R to better understand the numerical procedure. For those that prefer a non-R solution, the authors of WMM provide free software to calculate magnetic field, which can be found [here](https://www.ngdc.noaa.gov/geomag/WMM/soft.shtml).

# Usage

Currently, the only exported function is `GetMagneticFieldWMM`, which returns the orthogonal vector components of the main magnetic field (in nT) and secular variation field (in nT/yr) as predicted by WMM.

Example usage: 

1. Calculate expected magnetic field components at a benchmark location, mid 2017. Using the default value for WMM version, this will use the more recent, "out of cycle" coefficients from WMM2015v2 as opposed to the older WMM2015.
``` r

GetMagneticFieldWMM(
  lon = 240,
  lat = -80,
  height = 1e5,
  time = 2017.5
)
# $x
# [1] 5674.898
# 
# $y
# [1] 14793.08
# 
# $z
# [1] -50179.48
# 
# $xDot
# [1] 25.0788
# 
# $yDot
# [1] 1.523768
# 
# $zDot
# [1] 82.51838
```

2. Repeat the last calculation but use the older coefficients from WMM2015 that were replaced by WMM2015v2.
``` r
GetMagneticFieldWMM(
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
# $zDots
# [1] 86.24356
```

# Citations

1. Chulliat, A., W. Brown, P. Alken, S. Macmillan, M. Nair, C. Beggan, A. Woods,
B. Hamilton, B. Meyer and R. Redmon, 2019, Out-of-Cycle Update of the
US/UK World Magnetic Model for 2015-2020: Technical Note, National
Centers for Environmental Information, NOAA. doi: 10.25921/xhr3-0t19

2. Chulliat, A., S. Macmillan, P. Alken, C. Beggan, M. Nair, B. Hamilton, A.
Woods, V. Ridley, S. Maus and A. Thomson, 2015, The US/UK World
Magnetic Model for 2015-2020: Technical Report, National Geophysical
Data Center, NOAA. doi: 10.7289/V5TB14V7

3. Maus, S., S. Macmillan, S. McLean, B. Hamilton, A. Thomson,
M. Nair, and C. Rollins, 2010, The US/UK World Magnetic Model
for 2010-2015, NOAA Technical Report NESDIS/NGDC.

4. McLean, S., S. Macmillan, S. Maus, V. Lesur, A.
Thomson, and D. Dater, December 2004, The
US/UK World Magnetic Model for 2005-2010,
NOAA Technical Report NESDIS/NGDC-1. 

5. Macmillian, S. and J. M. Quinn, 2000. 
“The Derivation of the World Magnetic Model 2000,” 
British Geological Survey Technical Report WM/00/17R.

# Project Layout

The project has the following structure:

```
wmm/
├── .gitignore
├── .Rbuildignore
├── .Rprofile
├── .travis.yml
├── DESCRIPTION
├── LICENSE
├── NAMESPACE
├── README.md
├── wmm.Rproj
│   
├─── data-raw/
│   |── creatingSysData.R
│       
├───inst/
│   └── extdata/
│       |── WMM2000.csv
│       |── WMM2005.csv
│       |── WMM2010.csv
│       |── WMM2015.csv
│       |── WMM2015v2.csv
│       |── WMMTestValues.csv
│           
├── man/
├── R/
│   |── functions_coefficients.R
│   |── functions_geocoordinates.R
│   |── functions_legendre.R
│   |── functions_misc.R
│   |── functions_wmm.R
│   |── sysdata.rda
│       
└── tests/
    │── testthat.R
    │   
    └── testthat/
        |── testWMM.R
```
