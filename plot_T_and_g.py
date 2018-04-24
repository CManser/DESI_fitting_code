import numpy as np
import matplotlib.pyplot as plt
model = 'pier'
path = '/Users/christophermanser/Storage/PhD_files/DESI/2percent/'
f = path + 'WD_list_'+model+'.lst'
x = np.arange(5000, 100000, 1000)
name, T1, T1_err, g1, g1_err, T2, T2_err, g2, g2_err, rv, StN = np.genfromtxt(f, unpack = True)
desi_name, desi_T, desi_g, wd_type, mag_desi = np.genfromtxt(path + 'desi_WD_info.dat', unpack = True)
plt.plot(x, x, color = 'black', ls = '--')
for i in range(name.size):
  print(name[i], desi_name[i])
  plt.scatter(desi_T[i][wd_type[i] == 1], T1[i][wd_type[i] == 1], color = 'blue', s=1, zorder = 100)
  plt.scatter(desi_T[i][wd_type[i] == 1], T2[i][wd_type[i] == 1], color = 'green', s=1)
  #plt.scatter(desi_T[wd_type == 0], T1[wd_type == 0], color = 'red')
plt.show()
plt.close()
for i in range(name.size):
  plt.scatter(desi_T[i][wd_type[i] == 1], 1 - (T1[i][wd_type[i] == 1])/(desi_T[i][wd_type[i] == 1]), color = 'blue', s=5, zorder = 100)
  plt.scatter(desi_T[i][wd_type[i] == 1], 1 - (T2[i][wd_type[i] == 1])/(desi_T[i][wd_type[i] == 1]), color = 'red', s=10, alpha = 0.1)
  plt.axhline(0, color = '0.5', ls = '--')
plt.show()
plt.close()


for i in range(name.size):
  print(name[i], desi_name[i])
  plt.scatter(desi_g[i][wd_type[i] == 1], g1[i][wd_type[i] == 1], color = 'blue', s=1, zorder = 100)
  plt.scatter(desi_g[i][wd_type[i] == 1], g2[i][wd_type[i] == 1], color = 'green', s=1)
x = np.arange(6, 11, 1)
plt.plot(x, x, color = 'black', ls = '--')
#plt.errorbar(g1[wd_type == 1], desi_g[wd_type == 1], g1_err[wd_type == 1], fmt = None)
#plt.scatter(desi_g[wd_type == 0], g1[wd_type == 0], color = 'red')
plt.show()
plt.close()


plt.scatter(desi_T, desi_g, s = 1, color = 'black')
plt.scatter(T1, g1, s = 1, color = 'red')
plt.show()
plt.close()