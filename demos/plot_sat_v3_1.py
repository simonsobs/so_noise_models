from __future__ import print_function

from so_models_v3 import SO_Noise_Calculator_Public_v3_1_1 as so_models

import matplotlib
matplotlib.use('pdf')
matplotlib.rc('font', family='serif', serif='cm10')
matplotlib.rc('text', usetex=True)
fontProperties = {'family':'sans-serif',
                  'weight' : 'normal', 'size' : 16}
import matplotlib.pyplot as plt

####################################################################
##                   demonstration of the code
####################################################################

mode=2
suffix='pdf'
fsky=0.1
ellmax=1e3

dset_label = 'SAT\\_V3.1'
survey = so_models.SOSatV3point1(mode)
corr_pairs = [(0,1),(2,3),(4,5)]

ylims = (1e-6,1e-2)
colors=['b','r','g','m','k','y']
corr_colors = ['orange', 'fuchsia', 'springgreen']

print(dset_label)
bands = survey.get_bands()
print("band centers: ", survey.get_bands(), "[GHz]")
print("beam sizes: "  , survey.get_beams(), "[arcmin]")
N_bands = len(bands)

#ell, N_ell_T_full,N_ell_P_full = survey.get_noise_curves(
ell, N_ell_T_full, N_ell_P_full = survey.get_noise_curves(
    fsky, ellmax, 1, full_covar=False, deconv_beam=True)

WN_levels = survey.get_white_noise(fsky)**.5

N_ell_P  = N_ell_P_full

print("white noise levels: "  , WN_levels, "[uK-arcmin]")

target = str(survey.__class__.__name__).split('.')[-1]

## plot the polarization noise curves
plt.figure()
for i in range(N_bands):
    plt.loglog(ell,N_ell_P[i], label='%i GHz (%s)' % (bands[i], dset_label),
               color=colors[i], ls='-', lw=2.)

# include correlated atmospheric noise across frequencies
for _c,(i,j) in []:#enumerate(corr_pairs):
    plt.loglog(ell, N_ell_P_full[i,j],
               label=r'$%i \times %i$ GHz atm.' % (bands[i],bands[j]),
               color=corr_colors[_c], lw=1.5)

plt.title(r"$N(\ell$) Polarization", fontsize=18)
plt.ylabel(r"$N(\ell$) [$\mu$K${}^2$]", fontsize=16)
plt.xlabel(r"$\ell$", fontsize=16)
plt.ylim(*ylims)
plt.xlim(10,500)
plt.legend(loc='upper left', ncol=2, fontsize=9)
plt.grid()
plt.savefig('%s_mode%i_fsky%.2f_SAT_P.%s' % (target, mode, fsky, suffix))
####################################################################
####################################################################


