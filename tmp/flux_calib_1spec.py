import numpy as np
import sys
import matplotlib.pyplot as plt
from scipy import optimize
import fitting_scripts_1spec as fitting_scripts
import scipy.interpolate


def spec_bin(w, f, e, bin_width = 40):
  """Bins a spectrum with a resolution of bin_width"""
  wav_range   = max(w) - min(w)
  num_of_bins = int( np.ceil(wav_range / bin_width) )
  wb   = np.zeros(num_of_bins)
  fb   = np.zeros(num_of_bins)
  eb   = np.zeros(num_of_bins)
  temp = min(w) + (bin_width*0.5)
  for i in range(num_of_bins):
    if f[(w >= (temp-bin_width*0.5)) & (w < (temp+bin_width*0.5))].size != 0:
      flux_range    = f[(w >= (temp-bin_width*0.5)) & (w < (temp+bin_width*0.5))]
      err_range     = e[(w >= (temp-bin_width*0.5)) & (w < (temp+bin_width*0.5))]
      wb[i], fb[i]  = temp, np.sum(flux_range*err_range**-2)/np.sum(err_range**-2)
      eb[i]         = np.sqrt(1/np.sum(err_range**2))
      eb[i]         = ( np.sqrt( np.sum( (err_range)**2 ) ) ) / err_range.size
    temp          = temp + bin_width
  #
  return wb[((wb != 0) & (fb != 0) & (eb != 0))], fb[((wb != 0) & (fb != 0) & (eb != 0))], eb[((wb != 0) & (fb != 0) & (eb != 0))]

def zordersclfct(y, e, t):
  """Calculates the zeroth order scaling factor for a
     model best fit. y, e and t MUST have the same wavelength
     scale!"""
  invar2 = e**-2
  Syt = sum(y*t*invar2)
  Stt = sum(t*t*invar2)
  scale_factor = Syt/Stt
  return scale_factor

print(sys.argv[1])
spectra = np.loadtxt(sys.argv[1],usecols=(0,1,2),unpack=True).transpose()
spectra = spectra[np.isnan(spectra[:,1])==False & (spectra[:,0]>3500)]
spec1, spec2, spec3 = spectra[:,0].copy(), spectra[:,1].copy(), spectra[:,2].copy()
#flux_uncorr_b = spectra[:,3].copy()
# Running the Flux calibration
resp_b, model = fitting_scripts.flux_calib(spectra, sys.argv[1])


# Plotting
fig = plt.figure(figsize= (10,5))
axes1 = fig.add_axes([0, 0.5, 1 , 0.5])
axes2 = fig.add_axes([0, 0, 1, 0.45])
#axes3 = fig.add_axes([0, 0, 0.45, 0.20])


axes2.text(4500, 0.3, 'T=%d\nlogg=%.2f'%(model[0], model[1]/100))
inp = sys.argv[1]
n = inp.split('/')[-1]
StN = np.sum( (spec2[(spec1 >= 4500.0) & (spec1 <= 4750.0)] / spec3[(spec1 >= 4500.0) & (spec1 <= 4750.0)] )) / spec2[(spec1 >= 4500.0) & (spec1 <= 4750.0)].size
axes2.text(  4500, 0.1, 'StN = %.1f'%(StN), bbox={'edgecolor':'white', 'facecolor':'white', 'alpha':1, 'pad':2}, fontsize = 8)

#######################################
s2_c, s3_c = spec2/resp_b[0], spec3/resp_b[0]

pspec2 = spec2/np.mean(spec2)
ps2_c = s2_c/np.mean(s2_c)
#pflux_uncorr_b = flux_uncorr_b/np.mean(flux_uncorr_b)
ps3_c = s3_c/np.mean(s2_c)

axes1.plot(spec1, pspec2, color = 'black', lw = 1)
axes1.plot(spec1, ps2_c, color = 'blue', alpha = 0.8, lw = 1)
#axes1.plot(spec1, pflux_uncorr_b, color = 'red', lw = 2)
axes2.plot(spec1, resp_b[1], color = 'black')#/max(r_tmp))
axes2.plot(spec1, resp_b[0], color = 'red')#/max(resp))
#axes2.plot(spec1, spec2)

#axes3.plot(spec1, ps2_c / pflux_uncorr_b)
#tmp = spec_bin(spec1, ps2_c/pflux_uncorr_b, ps3_c / pflux_uncorr_b)
#axes3.plot(tmp[0], tmp[1])
axes2.xaxis.set_ticklabels([])
axes1.set_xlim([3500, 10500])
axes2.set_xlim([3500, 10500])
axes1.set_ylim([-2, 4.5])
#axes3.set_xlim([3500, 10500])
axes2.set_ylim([-0.1, 2])
#axes3.set_ylim([0.7, 1.3])
axes2.axhline(0.95, color = 'black', ls = '--', lw = 2, zorder = 100)
axes2.axhline(1.05, color = 'black', ls = '--', lw = 2, zorder = 100)

tell_start = [6800, 7580, 8050, 8860, 9280]
tell_end   = [7000, 7700, 8350, 9150, 9780]

for i in range(len(tell_start)):
  axes1.axvspan(tell_start[i], tell_end[i], facecolor='#2ca02c', alpha=0.5)
  axes2.axvspan(tell_start[i], tell_end[i], facecolor='#2ca02c', alpha=0.5)

#######################################

sav_dir = '/Users/christophermanser/Storage/PhD_files/DESI/BOSS_testing/spPlate_flux_text/'

#Blue Save
sav_dat = np.transpose(np.vstack((spec1, s2_c, s3_c)))
specname_b = sys.argv[1][:-4].split('/')[-1].split('_')[-1]
sav_name = 'fluxcorr_' + specname_b + '.dat'
np.savetxt(sav_dir + sav_name, sav_dat)

plt.savefig(sav_dir + sys.argv[1][:-4].split('/')[-1] + '.png', dpi = 300, bbox_inches = 'tight')
#plt.show()
plt.close()