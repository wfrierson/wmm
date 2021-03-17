# wmm 1.1.0
1. Updated Gauss coefficients for WMM2020.
2. Significantly improved numerical stability by using closed-form equation for associated Legendre polynomials, instead of recursion.
3. Improved speed by ~7x by using multidimensional arrays.
4. Per WMM2020, this version now displays a blackout zone warning when the horizontal intensity is below 6000 nT.
5. In addition to orthogonal magnetic field components, `GetMagneticFieldWMM` returns the magnetic field elements _h_, _f_, _i_, & _d_ as well as their secular variation.

# wmm 1.0.0
First release of wmm package. Contains 1 exported function, `GetMagneticFieldWMM`, which calculates magnetic field.
