import numpy as np
import sys
import fitting_scripts
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
# BLUE
spectra = np.loadtxt(sys.argv[1],usecols=(0,1,2),unpack=True).transpose()
spectra = spectra[np.isnan(spectra[:,1])==False & (spectra[:,0]>3500)]
#spectra[:,2]=spectra[:,2]**-2
spec1, spec2, spec3 = spectra[:,0].copy(), spectra[:,1].copy(), spectra[:,2].copy()

# Running the Flux calibration
resp_b, resp_r, resp_z, model = fitting_scripts.flux_calib(spectra, sys.argv[1])
s2_c, s3_c = spec2/resp_b[0], spec3/resp_b[0]

#RED
inp_r = sys.argv[1][:-20] + 'r' + sys.argv[1][-19:]
s_r = np.loadtxt(inp_r,usecols=(0,1,2),unpack=True).transpose()
s_r = s_r[np.isnan(s_r[:,1])==False & (s_r[:,0]>3500)]
s_r1, s_r2, s_r3 = s_r[:,0].copy(), s_r[:,1].copy(), s_r[:,2].copy()
sr2_c, sr3_c = s_r2/resp_r[0], s_r3/resp_r[0]

#ZED
inp_z = sys.argv[1][:-20] + 'z' + sys.argv[1][-19:]
s_z = np.loadtxt(inp_z,usecols=(0,1,2),unpack=True).transpose()
s_z = s_z[np.isnan(s_z[:,1])==False & (s_z[:,0]>3500)]
s_z1, s_z2, s_z3 = s_z[:,0].copy(), s_z[:,1].copy(), s_z[:,2].copy()
sz2_c, sz3_c = s_z2/resp_z[0], s_z3/resp_z[0]

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


f_dm_interp = scipy.interpolate.interp1d(w_dm, f_dm, bounds_error = False)(spec1)
f_dm_interp_r = scipy.interpolate.interp1d(w_dm, f_dm, bounds_error = False)(s_r1)
f_dm_interp_z = scipy.interpolate.interp1d(w_dm, f_dm, bounds_error = False)(s_z1)


sav_dir = '/Users/christophermanser/Storage/PhD_files/DESI/flux_calibration_testing/calib_files/'

#Blue Save
sav_dat = np.transpose(np.vstack((spec1, spec2, spec3, resp_b[1], resp_b[0], s2_c, s3_c, f_dm_interp)))
specname_b = sys.argv[1][:-4].split('/')[-1][13:]
sav_name = 'final_' + specname_b + '.dat'
np.savetxt(sav_dir + sav_name, sav_dat)

#Red Save
sav_dat = np.transpose(np.vstack((s_r1, s_r2, s_r3, resp_r[1], resp_r[0], sr2_c, sr3_c, f_dm_interp_r)))
specname_r = specname_b.replace('-b', '-r')
sav_name = 'final_' + specname_r + '.dat'
np.savetxt(sav_dir + sav_name, sav_dat)

#Red Save
sav_dat = np.transpose(np.vstack((s_z1, s_z2, s_z3, resp_z[1], resp_z[0], sz2_c, sz3_c, f_dm_interp_z)))
specname_z = specname_b.replace('-b', '-z')
sav_name = 'final_' + specname_z + '.dat'
np.savetxt(sav_dir + sav_name, sav_dat)

