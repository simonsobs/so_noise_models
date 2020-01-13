from so_models_v3 import SO_Noise_Calculator_Public_v3_1_1 as so_models

import matplotlib
matplotlib.use('pdf')
matplotlib.rc('font', family='serif', serif='cm10')
matplotlib.rc('text', usetex=True)
fontProperties = {'family':'sans-serif',
                  'weight' : 'normal', 'size' : 16}
import matplotlib.pyplot as plt

####################################################################
####################################################################
##                   demonstration of the code
####################################################################

mode=2
suffix='pdf'
fsky=0.4
ellmax=1e4
el=50.

dset_label = 'LAT\\_V3.1'
lat = so_models.SOLatV3point1(mode, el=el)
corr_pairs = [(0,1),(2,3),(4,5)]

ylims = (5e-7,1e0)
colors=['b','r','g','m','k','y']
corr_colors = ['orange', 'fuchsia', 'springgreen']

print(dset_label)
bands = lat.get_bands()
print("band centers: ", lat.get_bands(), "[GHz]")
print("beam sizes: "  , lat.get_beams(), "[arcmin]")
N_bands = len(bands)

ell, N_ell_LA_T_full,N_ell_LA_P_full = lat.get_noise_curves(
    fsky, ellmax, 1, full_covar=True, deconv_beam=True)

WN_levels = lat.get_white_noise(fsky)**.5

N_ell_LA_T  = N_ell_LA_T_full[range(N_bands),range(N_bands)]
N_ell_LA_Tx = [N_ell_LA_T_full[i,j] for i,j in corr_pairs]
N_ell_LA_P  = N_ell_LA_P_full[range(N_bands),range(N_bands)]
N_ell_LA_Px = [N_ell_LA_P_full[i,j] for i,j in corr_pairs]

print("white noise levels: "  , WN_levels, "[uK-arcmin]")

target = str(lat.__class__.__name__).split('.')[-1]

## plot the temperature noise curves
plt.figure()
for i in range(N_bands):
    plt.loglog(ell,N_ell_LA_T[i], label='%i GHz (%s)' % (bands[i], dset_label),
               color=colors[i], ls='-', lw=2.)

# include correlated atmospheric noise across frequencies
for _c,(i,j) in enumerate(corr_pairs):
    plt.loglog(ell, N_ell_LA_T_full[i,j],
               label=r'$%i \times %i$ GHz atm.' % (bands[i],bands[j]),
               color=corr_colors[_c], lw=1.5)

plt.title(r"$N(\ell$) Temperature", fontsize=18)
plt.ylabel(r"$N(\ell$) [$\mu$K${}^2$]", fontsize=16)
plt.xlabel(r"$\ell$", fontsize=16)
plt.ylim(*ylims)
plt.xlim(100,10000)
plt.legend(loc='lower left', ncol=2, fontsize=8)
plt.grid()
plt.savefig('%s_mode%i_fsky%.2f_LAT_T.%s' % (target, mode, fsky, suffix))

## plot the polarization noise curves
plt.figure()
for i in range(N_bands):
    plt.loglog(ell,N_ell_LA_P[i], label='%i GHz (%s)' % (bands[i], dset_label),
               color=colors[i], ls='-', lw=2.)

# include correlated atmospheric noise across frequencies
for _c,(i,j) in enumerate(corr_pairs):
    plt.loglog(ell, N_ell_LA_P_full[i,j],
               label=r'$%i \times %i$ GHz atm.' % (bands[i],bands[j]),
               color=corr_colors[_c], lw=1.5)

plt.title(r"$N(\ell$) Polarization", fontsize=18)
plt.ylabel(r"$N(\ell$) [$\mu$K${}^2$]", fontsize=16)
plt.xlabel(r"$\ell$", fontsize=16)
plt.ylim(*ylims)
plt.xlim(100,10000)
plt.legend(loc='upper left', ncol=2, fontsize=9)
plt.grid()
plt.savefig('%s_mode%i_fsky%.2f_LAT_P.%s' % (target, mode, fsky, suffix))
####################################################################
####################################################################


