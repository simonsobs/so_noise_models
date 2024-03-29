SO lensing noise curves v3.1.1
==============================

This version containes the SO lensing noise curves in which the TB lensing noise curves are fixed compared to v3.1.0. 

The following document is the same as v3.1.0.

------------------------------

Generated by Toshiya Namikawa.
Based on v3.1.0 of the SO component separated noise curves.

These are minimum variance combinations of TT,TE,EE,EB and TB
with iterative reconstruction for fsky=0.4.


The columns are [ells, N_lensing_TT, N_lensing_TE, N_lensing_EE, N_lensing_TB, N_lensing_EB, N_lensing_Pol (EE+EB), N_lensing_MV (all), N_curl_TT, N_curl_TE, N_curl_EE, N_curl_TB, N_curl_EB, N_curl_Pol (EE+EB), N_curl_MV (all)]
"lensing" refers to the gradient mode -- this is usually what you need for science forecasts
"curl" refers to the curl mode -- this is used in null tests
TT, TE, EE, TB, EB are the various estimators. For temperature-only, use N_lensing_TT.
For polarization-only, use N_lensing_Pol (EE+EB).
For most science cases, you want to use N_lensing_MV (all).
All noise spectra are for lensing convergence; no factors of ell or 2pi.
i.e. these can be plotted directly against C_ell_kappa_kappa.

Files:
nlkk_v3_1_0{deproj0,deproj1,deproj2}_{SENS1,SENS2}_fsky0p4_it_lT30-3000_lP30-5000.dat

deproj0 -- no deprojection of any foregrounds (baseline)
deproj1 -- tSZ deprojected
deproj2 -- fiducial CIB SED deprojected

SENS1 -- baseline sensitivity (baseline)
SENS2 -- goal sensitivity


Baseline recommendation:
nlkk_v3_1_0deproj0_SENS1_fsky0p4_it_lT30-3000_lP30-5000.dat


