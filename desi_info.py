import numpy as np 
import astropy.io.fits as fits
import glob

file_path = '/Users/christophermanser/Storage/PhD_files/DESI/WD_ascii/'
path = '/Users/christophermanser/Storage/PhD_files/DESI/2percent/'
lst = glob.glob(file_path + '*b-arm.dat')
name = np.zeros(len(lst))
T_desi = np.zeros(len(lst))
g_desi = np.zeros(len(lst))
wd_type = np.zeros(len(lst))
mag_desi = np.zeros(len(lst))

hdulist = fits.open(path + 'truth.fits')
data = hdulist[1].data

desi_ID = data['TARGETID']
desi_T = data['TEFF']
desi_g = data['LOGG']
WD_type = data['TEMPLATESUBTYPE']
desi_mag = data['MAG']

for i in range(len(lst)):
  print(i, '/', len(lst))
  name[i] = lst[i][len(file_path):].split('_')[0]
  print(name[i])
  element = np.where(name[i] == desi_ID)
  T_desi[i] = desi_T[element]
  g_desi[i] = desi_g[element]
  mag_desi[i] = desi_mag[element]
  print(mag_desi[i])
  if WD_type[element][0] == 'DA': wd_type[i] = 1
  
sav_data = np.dstack((name, T_desi, g_desi, wd_type, mag_desi))[0]
np.savetxt(path + 'desi_wd_info.dat', sav_data)