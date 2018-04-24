import sys
import numpy as np
import matplotlib.pyplot as plt

print('\n' + sys.argv[1])
# Loads the input spectrum as sys.argv[1], first input, and normalises
spectra = np.loadtxt(sys.argv[1],usecols=(0,1,2),unpack=True).transpose()
spectra = spectra[np.isnan(spectra[:,1])==False & (spectra[:,0]>3500) & (spectra[:,2] != 0.0)]
spec_w = spectra[:,0]


StN = np.sum( spectra[:,1][(spec_w >= 4500.0) & (spec_w <= 4750.0)] / spectra[:,2][(spec_w >= 4500.0) & (spec_w <= 4750.0)] ) / spectra[:,1][(spec_w >= 4500.0) & (spec_w <= 4750.0)].size
print('StN = %.1f'%(StN))
x_1, y_1 = spectra[:,0][(spec_w >= 4500.0) & (spec_w <= 4750.0)], spectra[:,1][(spec_w >= 4500.0) & (spec_w <= 4750.0)]
coeff = np.polyfit(x_1, y_1, 3)
p = np.poly1d(coeff)(x_1)
y_2 = y_1 - p

y_2 = abs(y_2)
print(np.mean(y_1) / np.mean(y_2))

plt.plot(x_1, y_2)
plt.show()
plt.close()


plt.plot(spectra[:,0], spectra[:,1])
plt.show()
plt.close()