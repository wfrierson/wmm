<!-- badges: start -->
[![Travis build status](https://travis-ci.org/wfrierson/wmm.svg?branch=master)](https://travis-ci.org/wfrierson/wmm)
<!-- badges: end -->

# wmm
R implementation of World Magnetic Model (WMM). This is a work in progress.

__Note__: I am not associated with those who built any version of WMM. Instead, I have recreated their algorithm for others to use in R.

# Usage

Currently, the only exported function is GetMagneticFieldWMM, which returns the orthogonal vector components of the main magnetic field (in nT) as predicted by WMM.

Example usage: 

GetMagneticFieldWMM(
  lon = 240,
  lat = -80,
  height = 1e5,
  time = 2017.5,
  wmmVersion = 'WMM2015'
)

Expected output:

$X
[1] 5683.518

$Y
[1] 14808.85

$Z
[1] -50163.01

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
