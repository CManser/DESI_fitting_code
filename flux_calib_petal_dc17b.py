import numpy as np
import sys
import matplotlib.pyplot as plt
from scipy import optimize
import fitting_scripts
import scipy.interpolate

def zordersclfct(y, e, t):
  """Calculates the zeroth order scaling factor for a
     model best fit. y, e and t MUST have the same wavelength
     scale!"""
  invar2 = e**-2
  Syt = sum(y*t*invar2)
  Stt = sum(t*t*invar2)
  scale_factor = Syt/Stt
  return scale_factor

spectra = np.loadtxt(sys.argv[1],usecols=(0,1,2),unpack=True).transpose()
spectra = spectra[np.isnan(spectra[:,1])==False & (spectra[:,0]>3500)]
#spectra[:,2]=spectra[:,2]**-2
spec1, spec2, spec3 = spectra[:,0].copy(), spectra[:,1].copy(), spectra[:,2].copy()

# Running the Flux calibration
resp_b, resp_r, resp_z, model = fitting_scripts.flux_calib(spectra, sys.argv[1])





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

# Plotting
fig = plt.figure(figsize= (12,8))
axes1 = fig.add_axes([0.5, 0.69, 0.5 , 0.3])
axes2 = fig.add_axes([0  , 0.7633, 0.45, 0.2366])
axes3 = fig.add_axes([0, 0.69, 0.45, 0.0666])

axes4 = fig.add_axes([0.5, 0.34, 0.5 , 0.30])
axes5 = fig.add_axes([0  , 0.43, 0.45, 0.2366])
axes6 = fig.add_axes([0, 0.36333, 0.45, 0.0666])

axes7 = fig.add_axes([0.5, 0, 0.5 , 0.30])
axes8 = fig.add_axes([0  , 0.066666, 0.45, 0.2666])
axes9 = fig.add_axes([0, 0, 0.45, 0.0666])

axes2.text(4500, 0.5, 'T=%d\nlogg=%.2f'%(model[0], model[1]/100))
inp = sys.argv[1]
n = inp.split('/')[-1]
StN = np.sum( (spec2[(spec1 >= 4500.0) & (spec1 <= 4750.0)] / spec3[(spec1 >= 4500.0) & (spec1 <= 4750.0)] )) / spec2[(spec1 >= 4500.0) & (spec1 <= 4750.0)].size
axes2.text(  4500, 0.3, 'StN = %.1f'%(StN), bbox={'edgecolor':'white', 'facecolor':'white', 'alpha':1, 'pad':2}, fontsize = 8)

#######################################
s2_c, s3_c = spec2/resp_b[0], spec3/resp_b[0]
axes1.plot(spec1, spec2, color = 'black', lw = 1)
axes1.plot(spec1, s2_c, color = 'blue', alpha = 0.8, lw = 1)
f_dm_interp = scipy.interpolate.interp1d(w_dm, f_dm, bounds_error = False)(spec1)
f_dm_interp_scl = f_dm_interp * zordersclfct(s2_c, s3_c, f_dm_interp)
axes1.plot(spec1, f_dm_interp_scl, color = 'red', lw = 1)
axes2.plot(spec1, resp_b[1], color = 'black')#/max(r_tmp))
axes2.plot(spec1, resp_b[0], color = 'red')#/max(resp))
#axes2.plot(spec1, spec2)

axes3.plot(spec1, s2_c / f_dm_interp_scl)
axes2.xaxis.set_ticklabels([])
axes1.set_xlim([3500, 6000])
axes2.set_xlim([3500, 6000])
axes3.set_xlim([3500, 6000])
axes3.axhline(0.95, color = '0.5', ls = '--', lw = 0.5)
axes3.axhline(1.05, color = '0.5', ls = '--', lw = 0.5)
#######################################


#RED
inp_r = sys.argv[1][:-20] + 'r' + sys.argv[1][-19:]
s_r = np.loadtxt(inp_r,usecols=(0,1,2),unpack=True).transpose()
s_r = s_r[np.isnan(s_r[:,1])==False & (s_r[:,0]>3500)]
s_r1, s_r2, s_r3 = s_r[:,0].copy(), s_r[:,1].copy(), s_r[:,2].copy()
#######################################
sr2_c, sr3_c = s_r2/resp_r[0], s_r3/resp_r[0]
axes4.plot(s_r1, s_r2, color = 'black', lw = 1)
axes4.plot(s_r1, sr2_c, color = 'blue', alpha = 0.8, lw = 1)
f_dm_interp_r = scipy.interpolate.interp1d(w_dm, f_dm, bounds_error = False)(s_r1)
f_dm_interp_scl_r = f_dm_interp_r * zordersclfct(sr2_c, sr3_c, f_dm_interp_r)
axes4.plot(s_r1, f_dm_interp_scl_r, color = 'red', lw = 1)
axes5.plot(s_r1, resp_r[1], color = 'black')#/max(r_tmp))
axes5.plot(s_r1, resp_r[0], color = 'red')#/max(resp))
#axes2.plot(spec1, spec2)

axes6.plot(s_r1, sr2_c / f_dm_interp_scl_r)
axes5.xaxis.set_ticklabels([])
axes5.set_xlim([5600, 7800])
axes6.set_xlim([5600, 7800])
axes6.set_xlim([5600, 7800])
axes6.axhline(0.95, color = '0.5', ls = '--', lw = 0.5)
axes6.axhline(1.05, color = '0.5', ls = '--', lw = 0.5)
#######################################

#ZED
inp_z = sys.argv[1][:-20] + 'z' + sys.argv[1][-19:]
s_z = np.loadtxt(inp_z,usecols=(0,1,2),unpack=True).transpose()
s_z = s_z[np.isnan(s_z[:,1])==False & (s_z[:,0]>3500)]
s_z1, s_z2, s_z3 = s_z[:,0].copy(), s_z[:,1].copy(), s_z[:,2].copy()
#######################################
sz2_c, sz3_c = s_z2/resp_z[0], s_z3/resp_z[0]
axes7.plot(s_z1, s_z2, color = 'black', lw = 1)
axes7.plot(s_z1, sz2_c, color = 'blue', alpha = 0.8, lw = 1)
f_dm_interp_z = scipy.interpolate.interp1d(w_dm, f_dm, bounds_error = False)(s_z1)
f_dm_interp_scl_z = f_dm_interp_z * zordersclfct(sz2_c, sz3_c, f_dm_interp_z)
axes7.plot(s_z1, f_dm_interp_scl_z, color = 'red', lw = 1)
axes8.plot(s_z1, resp_z[1], color = 'black')#/max(r_tmp))
axes8.plot(s_z1, resp_z[0], color = 'red')#/max(resp))
#axes2.plot(spec1, spec2)

axes9.plot(s_z1, sz2_c / f_dm_interp_scl_z)
axes8.xaxis.set_ticklabels([])
axes7.set_xlim([7400, 9900])
axes8.set_xlim([7400, 9900])
axes9.set_xlim([7400, 9900])
axes9.axhline(0.95, color = '0.5', ls = '--', lw = 0.5)
axes9.axhline(1.05, color = '0.5', ls = '--', lw = 0.5)
#######################################





sav_dir = '/Users/christophermanser/Storage/PhD_files/DESI/flux_calibration_testing/flux_calib/'
plt.savefig(sav_dir + sys.argv[1][:-4].split('/')[-1] + '.png', dpi = 300, bbox_inches = 'tight')
#plt.show()
plt.close()