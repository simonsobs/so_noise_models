For questions, contact so_tac@simonsobservatory.org .

This directory contains effective noise power spectra for SO LAT component-separated CMB T, E, B, and Compton-y maps.



If you need T, y, E, or B noise curves and aren't too worried about foregrounds, use:

[ell] [N_ell^TT] [N_ell^yy]:
SOV3_atmv1_T_default_noisecurves_deproj0_SENS1_mask_16000_ell_TT_yy.txt [baseline SO noise]
SOV3_atmv1_T_default_noisecurves_deproj0_SENS2_mask_16000_ell_TT_yy.txt [goal SO noise]

[ell] [N_ell^EE] [N_ell^BB]:
SOV3_pol_default1-4-2_noisecurves_deproj0_SENS1_mask_16000_ell_EE_BB.txt [baseline SO noise]
SOV3_pol_default1-4-2_noisecurves_deproj0_SENS2_mask_16000_ell_EE_BB.txt [goal SO noise]




The B-mode noise here is only for the SO LAT, i.e., not for primordial gravitational wave science (SATs).

The methodology is described in Sec. 2 of the SO science forecasting paper: https://arxiv.org/abs/1808.07445

The T and y noise curves have been re-computed (since the paper was published) using the latest V3_1_0 noise model, which corrects a minor atmospheric noise bug.

The E- and B-mode noise curves are unchanged from those computed in the paper.

All curves computed here are for an observed sky fraction of fsky = 0.4 ("mask_16000" in the file names).

The other file name conventions are:
- "SENS1" = baseline SO noise levels
- "SENS2" = goal SO noise levels
- "deproj0" = standard ILC
- "deproj1" = constrained ILC
  * tSZ deprojected for CMB T
  * CMB deprojected for tSZ y
  * fiducial polarized dust SED deprojected for CMB E and B
- "deproj2" = constrained ILC
  * fiducial CIB SED deprojected for CMB T
  * fiducial CIB SED deprojected for tSZ y
  * fiducial polarized synchrotron deprojected for CMB E and B
- "deproj3" = doubly constrained ILC
  * tSZ and CIB deprojected for CMB T
  * CMB and CIB deprojected for tSZ y
  * fiducial polarized dust and synchrotron SEDs deprojected for CMB E and B

The CMB T, E, and B noise power spectra are in units of uK^2 (no factors of ell, ell+1, etc.).

The Compton-y noise power spectra are dimensionless.

The first column of all files is the multipole, ell.  The others are as labeled in the filenames.

For the inclusion of additional Planck information, as described in Sec. 2.6 of the SO science forecasting paper, we use the specifications given in Table IV of https://arxiv.org/abs/1509.07471 .
