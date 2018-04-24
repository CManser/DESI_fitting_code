import numpy as np 
import astropy.io.fits as fits

f = 'WD_list_correct_err.lst'
path = '/Users/christophermanser/Storage/PhD_files/DESI/2percent/'
name, T1, T1_err, g1, g1_err, T2, T2_err, g2, g2_err, rv = np.genfromtxt(f, unpack = True)
T_desi = np.zeros(name.size)
g_desi = np.zeros(name.size)
wd_type = np.zeros(name.size)
desi_mag = np.zeros(name.size)

hdulist = fits.open(path + 'truth.fits')
data = hdulist[1].data

desi_ID = data['TARGETID']
desi_T = data['TEFF']
desi_g = data['LOGG']
WD_type = data['TEMPLATESUBTYPE']
desi_mag = data['MAG']

for i in range(name.size):
  element = np.where(name[i] == desi_ID)
  T_desi[i] = desi_T[element]
  g_desi[i] = desi_g[element]
  print(T1[i], desi_T[element])
  print(g1[i], desi_g[element])
  print(WD_type[element][0])
  if WD_type[element][0] == 'DA': wd_type[i] = 1
  
sav_data = np.dstack((name, T1, T1_err, g1, g1_err, T2, T2_err, g2, g2_err, rv, T_desi, g_desi, wd_type, desi_mag))[0]
np.savetxt(path + 'fit_vs_desi_model_correct_err.dat', sav_data)