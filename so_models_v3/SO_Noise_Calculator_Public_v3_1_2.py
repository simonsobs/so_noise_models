"""Simons Observatory LAT Noise Model

This is release v3.1.2.


Version 3.1.0 brings an update to the LAT T noise model, relative to
v3.0.  Version 3.1.1 brings the SAT noise model into the same file and
framework, for convenience (but the SAT noise is the same as for 3.0).

This code includes one SO SAT noise model:

  SOSatV3point1

This code includes two SO LAT noise models:

  SOLatV3
  - the original V3 noise model, but expanded to include elevation
    dependence. This should reproduced original V3 results at the
    reference elevation of 50 degrees.

  SOLatV3point1
  - like SOLatV3 but with an updated atmospheric 1/f model.

This code is based on the original SO LAT Noise Model released with
the SO forecasts paper, but with some functional changes:

- Object-oriented organization to make it easier to swap in different
  noise models for the same calculation.
- Beam deconvolution is optional.
- Plotting code is not included, and instead found separately in this
  repo.

Version 3.1.2 fixes the noise elevation scaling for the LAT MF bands.
The MF noise elevation scalings were using an average of the feedhorn
and lenslet values, but we are only using feedhorns at MF now.
  
"""

from __future__ import print_function

import numpy as np

def get_atmosphere_params(freqs, version=1, el=None):
    """Returns atmospheric noise power parameters, for an ACTPol optics
    tube.

    Arguments:

      freqs: array of frequencies (in GHz) to process.  This function only
        handles the standard SO frequency values.
      version: version of the C factors to return.  See below.
      el: elevation angle, in degrees.  Default is 50.

    Returns (C_array, alpha_array), where each array has the same
    shape as freqs.  The red noise contribution to each band is then:

          N_red(ell) =   C * (ell/1000)^alpha

    with units of [uK^2 sec].

    The model is naturally calibrated for boresight elevation el = 50
    degrees.  A simple rescaling (proportional to csc(el)) is applied
    for other values of el.

    In the present model, alpha=-3.5 always.

    version=0: This atmospheric model was used in SO V3 forecasts but
    contains an error.

    version=1: This atmospheric model is better than version=0, in
    that at least one error has been corrected.  The relative
    atmospheric powers have been measured in AT model, and then
    calibrated to ACT.  Low frequency results are inflated somewhat
    because ACT sees more power at 90 GHz than predicted by this
    modeling.
    """
    # This el_correction is quite naive; atmosphere may be less
    # structured (i.e. smoothed out) at lower elevation because of the
    # increased thickness.
    if el is None:
        el = 50.
    el_correction = np.sin(50.*np.pi/180) / np.sin(el*np.pi/180)
    if version == 0:
        data_bands = np.array([ 27., 39., 93., 145., 225., 280.])
        data_C = np.array([200., 77., 1800., 12000., 68000., 124000.])
        data = {}
        for b,C in zip(data_bands, data_C):
            data[b] = C
        # Jam in a fake value for 20 GHz.
        data[20.] = 200.
    elif version == 1:
        data_bands = np.array([ 20., 27., 39., 93., 145., 225., 280.])
        data_C = np.array([  3.45565594e+02 * 1.4,
                             6.00734484e+01 * 1.4,
                             1.45291936e+01 * 1.4,
                             6.81426710e+02 * 1.4,
                             1.20000000e+04,
                             2.66015359e+05,
                             1.56040514e+06,])
        data = {}
        for b,C in zip(data_bands, data_C):
            data[b] = C
    return (np.array([data[f] * el_correction**2 for f in freqs]),
            np.array([-3.5 for f in freqs]))

def rolloff(ell, ell_off=None, alpha=-4, patience=2.):
    """Get a transfer function T(ell) to roll off red noise at ell <
    ell_off.  ell should be an ndarray.  Above the cut-off,
    T(ell>=ell_off) = 1.  For T(ell<ell_off) will roll off smoothly,
    approaching T(ell) \propto ell^-alpha.  The speed at which the
    transition to the full power-law occurs is set by "patience";
    patience (> 1) is the maximum allowed value of:

                       T(ell) * ell**alpha
                 -----------------------------
                  T(ell_off) * ell_off**alpha

    I.e, supposing you were fighting an ell**alpha spectrum, the
    roll-off in T(ell) will be applied aggressively enough that
    T(ell)*ell**alpha does not rise above "patience" times its value
    at ell_off.

    """
    if ell_off is None or ell_off <= 0:
        return np.ones(ell.shape)
    L2 = ell_off
    L1 = L2 * patience ** (2./alpha)
    x = -np.log(ell / L2) / np.log(L1 / L2)
    beta = alpha * np.log(L1 / L2)
    output = x*0
    output[x<0]  = (-x*x)[x<0]
    output[x<-1] = (1 + 2*x)[x<-1]
    return np.exp(output * beta)


"""See lower for subclasses of SOTel -- instrument specific
parameters are configured in __init__.  """


class SOTel:
    """Base class for SO SAT and LAT noise models.  The sub-class
    constructors should set up all the internal variables and then
    call precompute().  Then the noise levels can be obtained with
    get_white_noise / get_noise_curves.

    """
    # Switches that tell get_noise_curves what to return.
    has_T = False
    has_P = False

    # In-tube, cross-band correlation coefficients for red ("atmos").
    Tatmos_band_corr = 0.
    Patmos_band_corr = 0.

    # Factor by which to attenuate LAT atmospheric power, given FOV
    # relative to ACT?
    Tatmos_FOV_mod = 0.5

    def __init__(self, *args, **kwargs):
        raise RuntimeError('You should subclass this.')

    @property
    def n_bands(self):
        return len(self.bands) if self.bands is not None else 0

    def get_bands(self):
        return self.bands.copy()

    def get_beams(self):
        return self.beams.copy()

    def precompute(self, N_tubes, N_tels=1):

        white_noise_el_rescale = np.array([1.] * len(self.bands))
        if self.el is not None:
            el_data = self.el_noise_params
            el_lims = el_data.get('valid')
            if el_lims[0] == 'only':
                assert(self.el == el_lims[1])  # noise model only valid at one elevation...
            else:
                assert(el_lims[0] <= self.el) and (self.el <= el_lims[1])
                band_idx = np.array([np.argmin(abs(el_data['bands'] - b)) for b in self.bands])
                assert(np.all(abs(np.array(el_data['bands'])[band_idx] - self.bands) < 5))
                coeffs = el_data['coeffs']
                white_noise_el_rescale = np.array(
                    [el_noise_func(coeffs[i], self.el) / el_noise_func(coeffs[i], 50.)
                     for i in band_idx])

        # Accumulate total white noise level and atmospheric
        # covariance matrix for this configuration.
        
        band_weights = np.zeros(self.n_bands)
        for (tube_name, tube_count) in N_tubes:
            tube_noise = self.tube_configs[tube_name] * white_noise_el_rescale
            s = (tube_noise != 0)
            band_weights[s] += tube_count * N_tels * tube_noise[s]**-2

        self.band_sens = np.zeros(self.n_bands) + 1e9
        s = (band_weights > 0)
        self.band_sens[s] = band_weights[s]**-0.5

        # Special for atmospheric noise model.
        C, alpha = get_atmosphere_params(self.bands, self.atm_version, el=self.el)
        self.Tatmos_C = C * self.Tatmos_FOV_mod
        self.Tatmos_alpha = alpha
        self.Tatmos_ell = 1000. + np.zeros(self.n_bands)

        # Compute covariant weight matrix (atmosphere parameters).
        cov_weight = np.zeros((self.n_bands,self.n_bands))
        pcov_weight = np.zeros((self.n_bands,self.n_bands))
        for (tube_name, tube_count) in N_tubes:
            # Get the list of coupled bands; e.g. [1,2] for MF.
            nonz = self.tube_configs[tube_name].nonzero()[0]
            for i in nonz:
                for j in nonz:
                    assert(cov_weight[i,j] == 0.) # Can't do overlapping
                                                  # tubes without weights.
                    T_corr = {True: 1., False: self.Tatmos_band_corr}[i==j]
                    cov_weight[i,j] += ( tube_count * N_tels /
                                         (T_corr *
                                          (self.Tatmos_C[i] * self.Tatmos_C[j])**.5) )
                    P_corr = {True: 1., False: self.Patmos_band_corr}[i==j]
                    pcov_weight[i,j] = P_corr

        # Reciprocate non-zero elements.
        s = (cov_weight!=0)
        self.Tatmos_cov = np.diag([1e9]*self.n_bands)
        self.Tatmos_cov[s] = 1./cov_weight[s]

        # Polarization is simpler...
        self.Patmos_cov = pcov_weight

    def get_survey_time(self):
        """Returns the effective survey time (survey_years * efficiency), in
        seconds.

        """
        t = self.survey_years * 365.25 * 86400.    ## convert years to seconds
        return t * self.survey_efficiency

    def get_survey_spread(self, f_sky, units='arcmin2'):
        """Returns the dilution factor that converts array instrument
        sensitivity (in units uK^2 sec) to map white noise level
        (units uK^2 arcmin^2).  Units are arcmin^2 / second (unless
        units='sr' is set)

        """
        # Factor that converts uK^2 sec -> uK^2 arcmin^2.
        A = f_sky * 4*np.pi
        if units == 'arcmin2':
            A *= (60*180/np.pi)**2
        elif units != 'sr':
            raise ValueError("Unknown units '%s'." % units)
        return A / self.get_survey_time()

    def get_white_noise(self, f_sky, units='arcmin2'):
        """Returns the survey white noise level, in temperature, for each
        band, in uK^2 arcmin2, for the specified f_sky (0 < f_sky <= 1).

        Pass units='sr' to get uK^2 steradian units.

        """
        return self.band_sens**2 * self.get_survey_spread(f_sky, units=units)

    def get_noise_curves(self, f_sky, ell_max, delta_ell, deconv_beam=True,
                         full_covar=False, rolloff_ell=None):
        """Get the noise curves N(ell) for all bands.

        The ell vector is determined by ell_max and delta_ell: ell =
        range(2, ell_max, delta_ell).

        The f_sky is the area of the survey in units of a full sky; (0
        < f_sky <= 1).

        Returns (ell, T_noise, P_noise).  If a model does not describe
        one of these spectra (has_T == False, or has_P == False), the
        corresponding spectrum will return as None.  Otherwise, the
        shape of T_noise and P_noise will be (n_bands, n_ell) if
        full_covar is False, and (n_bands, n_bands, n_ell) if
        full_covar is True.

        If deconv_beam is True, then the beam transfer functions are
        deconvolved, to give the effective noise level relative to a
        signal at each ell.

        If rolloff_ell is specified, a transfer function is applied to
        reduce red noise below this cutoff.  The transfer function at
        ell > rolloff_ell will be 1.  See code if you care about what
        happens below that.

        """
        ell = np.arange(2, ell_max, delta_ell)
        W = self.band_sens**2

        # Get full covariance; indices are [band,band,ell]
        ellf = (ell/self.Tatmos_ell[:,None])**(self.Tatmos_alpha[:,None])
        T_noise = self.Tatmos_cov[:,:,None] * (ellf[:,None,:] * ellf[None,:,:])**.5

        # P noise is tied directly to the white noise level.
        P_low_noise = (2*W[:,None]) * (ell / self.Patmos_ell[:,None])**self.Patmos_alpha[:,None]
        P_noise = (self.Patmos_cov[:,:,None] *
                   (P_low_noise[:,None,:] * P_low_noise[None,:,:])**.5)

        # Add in white noise on the diagonal.
        for i in range(len(W)):
            T_noise[i,i] += W[i]
            P_noise[i,i] += W[i] * 2

        if rolloff_ell is not None:
            # Use the same simple rolloff for all bands, T & P.
            gain = rolloff(ell, rolloff_ell)
            T_noise *= gain
            P_noise *= gain

        # Deconvolve beams.
        if deconv_beam:
            beam_sig_rad = self.get_beams() * np.pi/180/60 / (8.*np.log(2))**0.5
            beams = np.exp(-0.5 * ell*(ell+1) * beam_sig_rad[:,None]**2)
            T_noise /= (beams[:,None,:] * beams[None,:,:])
            P_noise /= (beams[:,None,:] * beams[None,:,:])

        # Diagonal only?
        if not full_covar:
            ii = range(self.n_bands)
            T_noise = T_noise[ii,ii]
            P_noise = P_noise[ii,ii]

        T_out, P_out = None, None
        if self.has_T:
            T_out = T_noise * self.get_survey_spread(f_sky, units='sr')
        if self.has_P:
            P_out = P_noise * self.get_survey_spread(f_sky, units='sr')

        return (ell, T_out, P_out)


""" Elevation-dependent noise model for LAT is captured below.  Detector
white noise is determined by a fixed, instrumental component as well
as a contribution from atmospheric loading, which is should scale
roughly with 1/sin(elevation).  The 'coeffs' below give coefficents A,
B for the model

   sens \propto A + B / sin(elevation)

For 27 GHz and higher, results are from the Simons Observatory LAT
noise model provided by Carlos Sierra and Jeff McMahon on March 11,
2019.  For 20 GHz, fits from Ben Racine and Denis Barkats are used
(for a slightly different instrument configuration that nevertheless
gives compatible results across all bands).  """

def el_noise_func(P, el):
    a, b = P
    return a + b / np.sin(el*np.pi/180)

SO_el_noise_func_params = {
    'threshold': {
        'valid': ('only', 50.),
    },
    'baseline': {
        'valid': (25., 70.),
        'bands': [20,27,39,93,145,225,280],
        'coeffs': [
            (178.59719595, 33.72945249),  # From B. Racine & D. Barkats.
            (.87, .09),                   # From Carlos Sierra and J. McMahon vvv
            (.64, .25),
            (.80, .14),
            (.76, .17),
            (.58, .30),
            (.49, .36),
        ],
    },
    'goal': {
        'valid': (25., 70.),
        'bands': [20,27,39,93,145,225,280],
        'coeffs': [
            (178.59719595, 33.72945249),  # From B. Racine & D. Barkats.
            (.85, .11),                   # From Carlos Sierra and J. McMahon vvv
            (.65, .25),
            (.76, .17),
            (.73, .19),
            (.55, .32),
            (.47, .37),
        ],
    },
}


class SOLatV3(SOTel):
    has_T = True
    has_P = True
    atm_version = 0
    Tatmos_band_corr = 0.9
    Patmos_band_corr = 0.9

    def __init__(self, sensitivity_mode=None, N_tubes=None, survey_years=5.,
                 survey_efficiency = 0.2*0.85, el=None):
        """Arguments:

          sensitivity_mode (int or string): Should be 'threshold',
            'baseline', or 'goal'.  Alternately you can pass 0, 1, or
            2.

          N_tubes: A list of tuples giving the survey-averaged number
            of each LAT tube in operation.  For example, the default
            is [('LF', 1), ('MF', 4), ('UHF', 2)], populating a total
            of 7 tubes in this LAT.  Fractional tubes are acceptable
            (imagine a tube were swapped out part way through the
            survey).

          survey_years: Total calendar years that the survey operates.

          survey_efficiency: Fraction of calendar time that may be
            used to compute map depth.

          el: Elevation, in degrees.  This affects white noise and red
            noise, through separate scalings.

        """
        # Define the instrument.
        self.bands = np.array([
            27., 39., 93., 145., 225., 280.])
        self.beams = np.array([
            7.4, 5.1, 2.2, 1.4, 1.0, 0.9])

        # Set defaults for survey area, time, efficiency
        self.survey_years = survey_years
        self.survey_efficiency  = survey_efficiency
        
        # Translate integer to string mode; check it.
        sens_modes = {0: 'threshold',
                      1: 'baseline',
                      2: 'goal'}
        if sensitivity_mode is None:
            sensitivity_mode = 'baseline'
        elif sensitivity_mode in sens_modes.keys():
            sensitivity_mode = sens_modes.get(sensitivity_mode)

        assert(sensitivity_mode in sens_modes.values())  # threshold,baseline,goal? 0,1,2?
        self.sensitivity_mode = sensitivity_mode
        
        # Sensitivities of each kind of optics tube, in uK rtsec, by
        # band.  0 represents 0 weight, not 0 noise!  Note the weird
        # scalings for MF (4**.5) and UHF (2**.5) are to convert to
        # per-tube values from the reference design sensitivity values
        # that included 4 MF and 2 UHF tubes.
        nar = np.array
        self.tube_configs = {
            'threshold': {
                'LF':  nar([  61.,  30.,    0,    0,    0,    0 ]),
                'MF':  nar([    0,    0,  6.5,  8.1,    0,    0 ])*4**.5,
                'UHF': nar([    0,    0,    0,    0,  17.,  42. ])*2**.5,
            },
            'baseline': {
                'LF':  nar([  48.,  24.,    0,    0,    0,    0 ]),
                'MF':  nar([    0,    0,  5.4,  6.7,    0,    0 ])*4**.5,
                'UHF': nar([    0,    0,    0,    0,  15.,  36. ])*2**.5,
            },
            'goal': {
                'LF':  nar([  35.,  18.,    0,    0,    0,    0 ]),
                'MF':  nar([    0,    0,  3.9,  4.2,    0,    0 ])*4**.5,
                'UHF': nar([    0,    0,    0,    0,  10.,  25. ])*2**.5,
            },
        }[sensitivity_mode]

        # The reference tube config.
        ref_tubes = [('LF', 1), ('MF', 4), ('UHF', 2)]

        if N_tubes is None:
            N_tubes = ref_tubes
        else:
            N_tubes = [(b,x) for (b,n),x in zip(ref_tubes, N_tubes)]
            
        ##
        ## T noise
        ##

        # Elevation stuff
        self.el_noise_params = SO_el_noise_func_params[sensitivity_mode]
        self.el = el

        ##
        ## P noise
        ##
        # "atmos" just means low-ell here.  We do not claim it
        # originates from atmospheric sources.
        self.Patmos_ell = 700. + np.zeros(self.n_bands)
        self.Patmos_alpha = -1.4 + np.zeros(self.n_bands)

        # Do general computations.
        self.precompute(N_tubes)


class SOLatV3point1(SOLatV3):
    atm_version = 1


# SAT support

class SOSatV3point1(SOTel):
    has_T = False
    has_P = True
    atm_version = 1

    def __init__(self, sensitivity_mode=None, N_tubes=None,
                 survey_years=5.,
                 survey_efficiency = 0.2*0.85, el=None,
                 one_over_f_mode=0):
        """Arguments:

          sensitivity_mode (int or string): Should be 'threshold',
            'baseline', or 'goal'.  Alternately you can pass 0, 1, or
            2.

          N_tubes: A list of tuples giving the survey-averaged number
            of each SAT type in operation.  For example, the default
            is [('LF', .4), ('MF', 1.6), ('UHF', 1)], which can be
            interpreted as 3 total instruments; 1 UHF instrument, 1 MF
            instrument, and one instrument that spends 60% of the
            survey as an MF and 40% of the survey as an LF."

          survey_years: Total calendar years that the survey operates.

          survey_efficiency: Fraction of calendar time that may be
            used to compute map depth.

          el: Elevation, in degrees.  The present SAT model does not
            support this parameter.

          one_over_f_mode: 0 or 1 to select 'pessimistic' or
            'optimistic' red-noise behavior, respectively.
        """
        # Define the instrument.
        self.bands = np.array([
            27., 39., 93., 145., 225., 280.])
        self.beams = np.array([
            91., 63., 30., 17., 11., 9.])

        # Set defaults for survey area, time, efficiency
        self.survey_years = survey_years
        self.survey_efficiency  = survey_efficiency
        
        # Translate integer to string mode; check it.
        sens_modes = {0: 'threshold',
                      1: 'baseline',
                      2: 'goal'}
        if sensitivity_mode is None:
            sensitivity_mode = 'baseline'
        elif sensitivity_mode in sens_modes.keys():
            sensitivity_mode = sens_modes.get(sensitivity_mode)

        assert(sensitivity_mode in sens_modes.values())  # threshold,baseline,goal? 0,1,2?
        self.sensitivity_mode = sensitivity_mode
        
        # Sensitivities of each kind of optics tube, in uK rtsec, by
        # band.  0 represents 0 weight, not 0 noise!  Note the weird
        # scaling for MF (2**.5) is to convert to per-tube values from
        # the reference design sensitivity values that included 2 MF
        # instruments.
        nar = np.array
        self.tube_configs = {
            'threshold': {
                'LF':  nar([  32.,  17.,    0,    0,    0,    0 ]),
                'MF':  nar([    0,    0,  4.6,  5.5,    0,    0 ])*2**.5,
                'UHF': nar([    0,    0,    0,    0,  11.,  26. ]),
            },
            'baseline': {
                'LF':  nar([  21.,  13.,    0,    0,    0,    0 ]),
                'MF':  nar([    0,    0,  3.4,  4.3,    0,    0 ])*2**.5,
                'UHF': nar([    0,    0,    0,    0,  8.6,  22. ]),
            },
            'goal': {
                'LF':  nar([  15.,  10.,    0,    0,    0,    0 ]),
                'MF':  nar([    0,    0,  2.4,  2.7,    0,    0 ])*2**.5,
                'UHF': nar([    0,    0,    0,    0,  5.7,  14. ]),
            },
        }[sensitivity_mode]

        # Save the elevation request.
        assert(el is None)  # Sorry, no SAT elevation function!
        self.el = el

        # The reference tube config.
        ref_tubes = [('LF', .4), ('MF', 1.6), ('UHF', 1)]

        if N_tubes is None:
            N_tubes = ref_tubes
        else:
            N_tubes = [(b,x) for (b,n),x in zip(ref_tubes, N_tubes)]
            
        ##
        ## T noise
        ##
        # This model does not describe the T noise.

        ##
        ## P noise
        ##
        # "atmos" just means low-ell here.  We do not claim it
        # originates from atmospheric sources.
        self.Patmos_alpha = np.array([-2.4,-2.4,-2.5,-3,-3,-3])
        if one_over_f_mode == 0:
            self.Patmos_ell = np.array([30.,30,50,50,70,100])
        elif one_over_f_mode == 1:
            self.Patmos_ell = np.array([15.,15,25,25,35,40])
        else:
            raise ValueError('Invalid one_over_f_mode')

        # Do general computations.
        self.precompute(N_tubes)
