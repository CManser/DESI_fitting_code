import fitting_scripts
import numpy as np
import sys
print('\n' + sys.argv[1])

spectra = np.loadtxt(sys.argv[1],usecols=(0,1,2),unpack=True).transpose()
spectra = spectra[np.isnan(spectra[:,1])==False & (spectra[:,0]>3500)]
fitting_scripts.fit_DA(spectra, plot = True)#,sys.argv[1])