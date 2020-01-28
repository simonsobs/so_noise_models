===============
so_noise_models
===============

This repository hosts publicly released N(ell) noise curve projection
code for the `Simons Observatory <https://simonsobservatory.org>`__.
The intention is that the full history of noise models will be
available here, to supplement published projections and simulations.

This main home of the repository is:
https://github.com/simonsobs/so_noise_models.

The noise model was originally described and applied in the
publication:

  Simons Observatory Collaboration.  "The Simons Observatory: science
  goals and forecasts".  JCAP 1902 (2019) 056.  arxiv:1808.07455

The paper may be found at the following links:

- https://dx.doi.org/10.1088/1475-7516/2019/02/056
- https://arxiv.org/abs/1808.07445


This repository is organized as follows:

- ``so_models_v3/``: Python 3 package providing the noise model code
  used in publications and publicly released simulations.  This
  consists of several independent sub-modules, representing each
  version of the noise code.  The usage of the models can vary
  substantially from version to version -- please consult code in
  ``demos/`` for typical usage patterns.
- ``demos/``: Code that demonstrates usage of the noise models, such
  as by producing noise curve plots.
- ``LAT_comp_sep_noise/`` - Effective noise power spectra for SO LAT
  component-separated CMB T, E, B, and Compton-y maps.  See dedicated
  README within.
- ``LAT_lensing_noise/`` - Lensing noise curves from SO LAT
  component-separated CMB T, E, B maps.  See dedicated
  README files within.

The ``so_models_v3`` package is pure python and thus can be imported
from the root level of this repository.  But you might want to install
it into your Python environment.  Conda users can simply run::

  pip3 install .

If you're not using Conda and want to install the package to your
local user package folder, run::

  python3 setup.py install --user
