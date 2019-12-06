"""Simons Observatory LAT Noise Model

This is release v3.1.0.  It is an update to the LAT noise model only.
v3.0.4 is still the most current SAT noise model.

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
- Plotting code is confined to __main__ block, so not run on simple
  import.

"""

from __future__ import print_function

import numpy as np

def get_atmosphere_C(freqs, version=1, el=None):
    """
    Returns atmospheric noise power at ell=1000, for an ACTPol optics
    tube.  In units of [uK^2 sec].  This only works for a few special
    frequencies.

    Basic model assumes el=50.  A simple rescaling (proportional to
    csc(el)) is applied for other values of el.

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
    return np.array([data[f] * el_correction**2 for f in freqs])


"""See lower for subclasses of SOLatType -- instrument specific
parameters are configured in __init__.  """

class SOLatType:
    def __init__(self, *args, **kwargs):
        raise RuntimeError('You should subclass this.')

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
        self.Tatmos_C = get_atmosphere_C(self.bands, self.atm_version,
                                         el=self.el) * self.FOV_mod
        self.Tatmos_ell = 1000. + np.zeros(self.n_bands)
        self.Tatmos_alpha = -3.5 + np.zeros(self.n_bands)

        # Compute covariant weight matrix (atmosphere parameters).
        cov_weight = np.zeros((self.n_bands,self.n_bands))
        pcov_weight = np.zeros((self.n_bands,self.n_bands))
        atm_rho = 0.9
        for (tube_name, tube_count) in N_tubes:
            # Get the list of coupled bands; e.g. [1,2] for MF.
            nonz = self.tube_configs[tube_name].nonzero()[0]
            for i in nonz:
                for j in nonz:
                    w = {True: 1., False: atm_rho}[i==j]
                    assert(cov_weight[i,j] == 0.) # Can't do overlapping
                                                  # tubes without weights.
                    cov_weight[i,j] += tube_count * N_tels / ( w * (
                        self.Tatmos_C[i] * self.Tatmos_C[j])**.5 )
                    pcov_weight[i,j] = w

        # Reciprocate non-zero elements.
        s = (cov_weight!=0)
        self.Tatmos_cov = np.diag([1e9]*self.n_bands)
        self.Tatmos_cov[s] = 1./cov_weight[s]

        # Polarization is simpler...
        self.Patmos_ell = 700. + np.zeros(self.n_bands)
        self.Patmos_alpha = -1.4 + np.zeros(self.n_bands)
        
        self.Patmos_cov = pcov_weight

    def get_survey_time(self):
        t = self.survey_years * 365.25 * 86400.    ## convert years to seconds
        return t * self.survey_efficiency

    def get_survey_spread(self, f_sky, units='arcmin2'):
        # Factor that converts uK^2 sec -> uK^2 arcmin^2.
        A = f_sky * 4*np.pi
        if units == 'arcmin2':
            A *= (60*180/np.pi)**2
        elif units != 'sr':
            raise ValueError("Unknown units '%s'." % units)
        return A / self.get_survey_time()

    def get_white_noise(self, f_sky, units='arcmin2'):
        return self.band_sens**2 * self.get_survey_spread(f_sky, units=units)

    def get_noise_curves(self, f_sky, ell_max, delta_ell, deconv_beam=True,
                         full_covar=False):
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

        sr_per_arcmin2 = (np.pi/180/60)**2
        return (ell,
                T_noise * self.get_survey_spread(f_sky, units='sr'),
                P_noise * self.get_survey_spread(f_sky, units='sr'))


""" Elevation-dependent noise model is parametrizerd below.  Detector
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
            ((.80+.74)/2, (.14+.19)/2),
            ((.76+.73)/2, (.17+.19)/2),
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
            ((.76+.69)/2, (.17+.22)/2),
            ((.73+.69)/2, (.19+.22)/2),
            (.55, .32),
            (.47, .37),
        ],
    },
}


class SOLatV3(SOLatType):
    atm_version = 0
    def __init__(self, sensitivity_mode=None, N_tubes=None, survey_years=5.,
                 survey_efficiency = 0.2*0.85, el=None):
        # Define the instrument.
        self.n_bands = 6
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
        # band.  0 represents 0 weight, not 0 noise...
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

        self.el_noise_params = SO_el_noise_func_params[sensitivity_mode]
        
        # Save the elevation request.
        self.el = el
        
        # Factor by which to attenuate atmospheric power, given FOV
        # relative to ACT?
        self.FOV_mod = 0.5

        # The reference tube config.
        ref_tubes = [('LF', 1), ('MF', 4), ('UHF', 2)]

        if N_tubes is None:
            N_tubes = ref_tubes
        else:
            N_tubes = [(b,x) for (b,n),x in zip(ref_tubes, N_tubes)]
            
        self.precompute(N_tubes)

    
class SOLatV3point1(SOLatV3):
    atm_version = 1
