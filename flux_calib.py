import numpy as np
import sys
import matplotlib.pyplot as plt
from scipy import optimize
import fitting_scripts

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
spectra[:,2]=spectra[:,2]**-0.5
spec1, spec2, spec3 = spectra[:,0], spectra[:,1], spectra[:,2]

#######################################
# For making an 'psuedo-UNCALIBRATED' spectrum
x1, x2, y1, y2 = 4500, 5000, 0.80, 0.90
tmp = ( ((y2-y1)/(x2-x1))*(spec1 - x1) + y1 ) # Straight line
cut = spec1 < x1
x = 1.5*(spec1[cut] - x1)/(x1 - spec1.min())
tmp[cut] = tmp[cut]*np.exp(-x*x) # Gaussian drop off
cut = spec1 > x2
x = 1.5*(spec1[cut] - x2)/(x2 - spec1.max())
tmp[cut] = tmp[cut]*np.exp(-x*x) # Gaussian drop off
spec2_adjust, spec3_adjust = spec2*tmp, spec3*tmp
spectra = np.stack((spec1,spec2_adjust,spec3_adjust), axis=-1)
#########################################
print('\n' + sys.argv[1])
resp, model = fitting_scripts.flux_calib(spectra, sys.argv[1])


fig = plt.figure(figsize= (12,6))
axes1 = fig.add_axes([0.5, 0, 0.5 , 1])
axes2 = fig.add_axes([0  , 0.2, 0.45, 0.8])
axes3 = fig.add_axes([0, 0, 0.45, 0.2])
axes1.plot(spec1, spec2,color = 'black', lw = 1)
axes2.text(4500, 0.5, '%d\n%.2f'%(model[0], model[1]/100))

#######################################
s2_c, s3_c = spec2_adjust/resp, spec3_adjust/resp
axes1.plot(spec1, spec2_adjust, color = 'red', lw = 1)
axes2.plot(spec1, resp/max(resp))
axes2.plot(spec1, tmp/max(tmp))
axes3.plot(spec1, (resp/max(resp)) / (tmp/max(tmp)))
axes2.set_xlim([3500, 6000])
axes3.set_xlim([3500, 6000])
axes3.axhline(0.95, color = '0.5', ls = '--', lw = 0.5)
axes3.axhline(1.05, color = '0.5', ls = '--', lw = 0.5)
#######################################
#s2_c, s3_c = spec2/resp, spec3/resp


axes1.plot(spec1, s2_c*zordersclfct(spec2[(spec1>4000) & (spec1<5500)], 
           spec3[(spec1>4000) & (spec1<5500)], s2_c[(spec1>4000) & (spec1<5500)]), 
           color = 'blue', alpha = 0.8, lw = 1)
sav_dir = '/Users/christophermanser/Storage/PhD_files/DESI/WD_ascii/WD_fits/20171117/flux_calib_test/'
plt.savefig(sav_dir + sys.argv[1][:-4].split('/')[-1] + '.png', dpi = 300, bbox_inches = 'tight')
#plt.show()
plt.close()