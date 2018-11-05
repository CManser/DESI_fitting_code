import numpy as np
import sys
import matplotlib.pyplot as plt
from scipy import optimize
import fitting_scripts_BOSS as fitting_scripts
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
spectra = np.loadtxt(sys.argv[1],usecols=(0,1,2,3,4),unpack=True).transpose()
spectra = spectra[np.isnan(spectra[:,1])==False & (spectra[:,0]>3500)]
spectra = spectra[spectra[:,2] != 0.0]
spectra[:,2]=spectra[:,2]**-0.5
spec1, spec2, spec3 = spectra[:,0].copy(), spectra[:,1].copy(), spectra[:,2].copy()
flux_uncorr_b = spectra[:,3].copy()
# Running the Flux calibration
resp_b, resp_r, model = fitting_scripts.flux_calib(spectra, sys.argv[1])




'''
# Finding the desi info and getting the model
file_name = sys.argv[1].split('/')[-1]
print(file_name)
obs = file_name.split('-')[-1][:8]
petal = file_name.split('-')[1][1:]
fiber = file_name.split('_')[-1][:-4]
fiber_big = str(int(fiber) + 500*int(petal))
check = obs + ' ' + fiber_big

list_path = '/Users/christophermanser/Storage/PhD_files/DESI/flux_calibration_testing/uncalibrated/'
wd_fibers = open(list_path + 'wd_fibers_teff_logg.dat', 'r')
lines = wd_fibers.readlines()
for i in lines:
  if check in i:
    fiber_info = i.split(' ')
    break
temp, logg = fiber_info[3], fiber_info[4]
dir_models = '/Users/christophermanser/Storage/PhD_files/DESI/WDFitting/WDModels_Koester.da2014_npy/'
model_file = 'da%06d_%d_2.7.npy'%(float(temp), float(logg)*100)
dm = np.load(dir_models+model_file)
w_dm, f_dm = dm[:,0], dm[:,1]

w_bm, f_bm, e_bm = model[2][:,0], model[2][:,1], model[2][:,2]
'''
# Plotting
fig = plt.figure(figsize= (12,8))
axes1 = fig.add_axes([0.5, 0.69, 0.5 , 0.3])
axes2 = fig.add_axes([0  , 0.7633, 0.45, 0.2366])
axes3 = fig.add_axes([0, 0.69, 0.45, 0.0666])

axes4 = fig.add_axes([0.5, 0.34, 0.5 , 0.30])
axes5 = fig.add_axes([0  , 0.43, 0.45, 0.2366])
axes6 = fig.add_axes([0, 0.36333, 0.45, 0.0666])

axes2.text(4500, 0.5, 'T=%d\nlogg=%.2f'%(model[0], model[1]/100))
inp = sys.argv[1]
n = inp.split('/')[-1]
StN = np.sum( (spec2[(spec1 >= 4500.0) & (spec1 <= 4750.0)] / spec3[(spec1 >= 4500.0) & (spec1 <= 4750.0)] )) / spec2[(spec1 >= 4500.0) & (spec1 <= 4750.0)].size
axes2.text(  4500, 0.3, 'StN = %.1f'%(StN), bbox={'edgecolor':'white', 'facecolor':'white', 'alpha':1, 'pad':2}, fontsize = 8)

#######################################
s2_c, s3_c = spec2/resp_b[0], spec3/resp_b[0]

pspec2 = spec2/np.mean(spec2)
ps2_c = s2_c/np.mean(s2_c)
pflux_uncorr_b = flux_uncorr_b/np.mean(flux_uncorr_b)
ps3_c = s3_c/np.mean(s2_c)

axes1.plot(spec1, pspec2, color = 'black', lw = 1)
axes1.plot(spec1, ps2_c, color = 'blue', alpha = 0.8, lw = 1)
axes1.plot(spec1, pflux_uncorr_b, color = 'red', lw = 2)
axes2.plot(spec1, resp_b[1], color = 'black')#/max(r_tmp))
axes2.plot(spec1, resp_b[0], color = 'red')#/max(resp))
#axes2.plot(spec1, spec2)

axes3.plot(spec1, ps2_c / pflux_uncorr_b)
tmp = spec_bin(spec1, ps2_c/pflux_uncorr_b, ps3_c / pflux_uncorr_b)
axes3.plot(tmp[0], tmp[1])
axes2.xaxis.set_ticklabels([])
axes1.set_xlim([3500, 6000])
axes2.set_xlim([3500, 6000])
axes3.set_xlim([3500, 6000])
axes2.set_ylim([-0.1, 2])
axes3.set_ylim([0.7, 1.3])
axes3.axhline(0.95, color = 'black', ls = '--', lw = 2, zorder = 100)
axes3.axhline(1.05, color = 'black', ls = '--', lw = 2, zorder = 100)
#######################################


#RED
inp_r = sys.argv[1][:-15] + 'r' + sys.argv[1][-14:]
s_r = np.loadtxt(inp_r,usecols=(0,1,2,3,4),unpack=True).transpose()
s_r = s_r[np.isnan(s_r[:,1])==False & (s_r[:,0]>3500)]
s_r = s_r[(s_r[:,2] != 0.0)]
s_r1, s_r2, s_r3 = s_r[:,0].copy(), s_r[:,1].copy(), s_r[:,2].copy()
s_r3 = s_r3**-0.5
flux_uncorr_r = s_r[:,3].copy()
#######################################
sr2_c, sr3_c = s_r2/resp_r[0], s_r3/resp_r[0]

ps_r2  = s_r2/np.mean(s_r2)
psr2_c = sr2_c/np.mean(sr2_c)
pflux_uncorr_r = flux_uncorr_r/np.mean(flux_uncorr_r)
psr3_c = sr3_c/np.mean(sr2_c)

axes4.plot(s_r1, ps_r2, color = 'black', lw = 1)
axes4.plot(s_r1, psr2_c, color = 'blue', alpha = 0.8, lw = 1)
axes4.plot(s_r1, pflux_uncorr_r, color = 'red', lw = 2)
axes5.plot(s_r1, resp_r[1], color = 'black')#/max(r_tmp))
axes5.plot(s_r1, resp_r[0], color = 'red')#/max(resp))
#axes2.plot(spec1, spec2)


axes6.plot(s_r1, psr2_c / pflux_uncorr_r)
tmp = spec_bin(s_r1, psr2_c/pflux_uncorr_r, psr3_c / pflux_uncorr_r)
axes6.plot(tmp[0], tmp[1])
axes5.xaxis.set_ticklabels([])
axes5.set_xlim([5600, 10000])
axes6.set_xlim([5600, 10000])
axes6.set_xlim([5600, 10000])
axes5.set_ylim([-0.1, 2])
axes6.set_ylim([0.75, 1.25])
axes6.axhline(0.95, color = 'black', ls = '--', lw = 2, zorder = 100)
axes6.axhline(1.05, color = 'black', ls = '--', lw = 2, zorder = 100)
#######################################

sav_dir = '/Users/christophermanser/Storage/PhD_files/DESI/BOSS_testing/wd_obs/flux_calib_test/'

#Blue Save
sav_dat = np.transpose(np.vstack((spec1, s2_c, s3_c)))
specname_b = sys.argv[1][:-4].split('/')[-1].split('_')[-1]
sav_name = 'fluxcorr_' + specname_b + '.dat'
np.savetxt(sav_dir + sav_name, sav_dat)

#Red Save
sav_dat = np.transpose(np.vstack((s_r1, sr2_c, sr3_c)))
specname_r = specname_b.replace('-b', '-r')
sav_name = 'fluxcorr_' + specname_r + '.dat'
np.savetxt(sav_dir + sav_name, sav_dat)

#plt.savefig(sav_dir + sys.argv[1][:-4].split('/')[-1] + '.png', dpi = 300, bbox_inches = 'tight')
plt.show()
plt.close()