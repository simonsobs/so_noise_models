Release date: 12 December 2019

For questions, contact so_tac@simonsobservatory.org .

This directory contains effective noise power spectra for SO LAT
component-separated CMB T, E, B, and Compton-y maps.

The B-mode noise here is only for the SO LAT, i.e., not for primordial
gravitational wave science (SATs).

The methodology is described in Sec. 2 of the SO science forecasting
paper: https://arxiv.org/abs/1808.07445

The T and y noise curves have been re-computed (since the paper was
published) using the latest V3_1_0 noise model, which corrects a minor
atmospheric noise bug.

The E- and B-mode noise curves are unchanged from those computed in
the paper.

All curves computed here are for an observed sky fraction of
fsky = 0.4, as labeled in the file names.

The other file name conventions are:
- "baseline" = baseline SO noise levels
- "goal" = goal SO noise levels
- "CMB" = CMB(+kSZ) blackbody component noise in units of uK^2
- "tSZ" = tSZ component noise in dimensionless Compton-y units

The file headers contain more information on the contents.

There are no factors of ell, ell+1, etc. applied to the spectra.

The first column of all files is the multipole, ell.  The others
are as labeled in the filenames.

For the inclusion of additional Planck information, as described in
Sec. 2.6 of the SO science forecasting paper, we use the
specifications given in Table IV of https://arxiv.org/abs/1509.07471.
