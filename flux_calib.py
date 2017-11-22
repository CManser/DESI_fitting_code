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

model_c='da2014'
basedir='/Users/christophermanser/Storage/PhD_files/DESI/WDFitting'
c = 299792.458 # Speed of light in km/s
plot = True
spectra = np.loadtxt(sys.argv[1],usecols=(0,1,2),unpack=True).transpose()
spectra = spectra[np.isnan(spectra[:,1])==False & (spectra[:,0]>3500)]
spectra[:,2]=spectra[:,2]**-0.5
spec1, spec2, spec3 = spectra[:,0], spectra[:,1], spectra[:,2]
print('\n' + sys.argv[1])
resp, model = fitting_scripts.flux_calib(sys.argv[1])
s2_c, s3_c = spec2/resp, spec3/resp
#s2_c, s3_c = spec2_adjust/resp, spec3_adjust/resp
print(model[0], model[1]/100)
fig = plt.figure(figsize= (12,6))
axes1 = fig.add_axes([0.5, 0, 0.5 , 1])
axes2 = fig.add_axes([0  , 0, 0.45, 1])
#axes1.plot(spec1, spec2_adjust, color = 'red', lw = 1)
axes1.plot(spec1, spec2,color = 'black', lw = 1)
axes1.plot(spec1, s2_c*zordersclfct(spec2[(spec1 > 4000) & (spec1 < 9000)], spec3[(spec1 > 4000) & (spec1 < 9000)], s2_c[(spec1 > 4000) & (spec1 < 9000)]), color = 'blue', alpha = 0.8, lw = 1)
print(resp)
axes2.text(5500, 0.5, str(model[0]) + ' ' + str(model[1]/100))
#sav_dir = '/Users/christophermanser/Storage/PhD_files/DESI/julian_SDSS_plates/'
#plt.savefig(sav_dir + sys.argv[1][-43:-4] + '_mockcalib.eps', bbox_inches = 'tight')
plt.show()
plt.close()