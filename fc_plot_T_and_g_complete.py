import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
model = 'da2014'
path = '/Users/christophermanser/Storage/PhD_files/DESI/flux_calibration_testing/'
f = path + 'fluxcalib_dc17b_all.lst'
x = np.arange(5000, 100000, 1000)
cm = plt.cm.get_cmap('viridis')
data = np.genfromtxt(f, unpack = True, dtype = None)

name = data['f0']
T1 = data['f1']
T1_err = data['f2']
g1 = data['f3']
g1_err = data['f4']
T2 = data['f5']
T2_err = data['f6']
g2 = data['f7']
g2_err = data['f8']
rv = data['f9']
StN = data['f10']



data2 = np.genfromtxt(path + 'uncalibrated/wd_fibers_teff_logg.dat', unpack = True, dtype = None)

obs_num = data2['f0']
big_fiber = data2['f1']
wd_type_tmp = data2['f2']
desi_T_tmp = data2['f3']
desi_g_tmp = data2['f4']
mag_desi_tmp = data2['f5']
flux_g = data2['f6']
flux_r = data2['f7']
flux_z = data2['f8']


desi_T = np.zeros(name.size)
desi_g = np.zeros(name.size)
desi_mag = np.zeros(name.size)
wd_type = np.zeros(name.size)

for i in range(name.size):
  print(i)
  tmp = name[i].split('-')[-1]
  obs = tmp[:8]
  fiber = int(tmp[-8:-4])
  petal = int(name[i].split('-')[-2][1])
  b_fiber = fiber + petal*500.0
  for j in range(obs_num.size):
    if (obs_num[j] == int(obs)) and (big_fiber[j] == int(b_fiber)):
      desi_T[i] = desi_T_tmp[j]
      desi_g[i] = desi_g_tmp[j]
      desi_mag[i] = mag_desi_tmp[j]
      if wd_type_tmp[j] == 'DA':
        wd_type[i] = 1
      else:
        wd_type[i] = 0
      break
      
#####################################################
matplotlib.rc('font',**{'size' : 18, 'family':'serif',
'serif':['serif']})
matplotlib.rc('text', usetex=True)
#matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['ps.usedistiller'] = 'xpdf'


######################################################
##################################################
##############     Plotting      #################
##################################################
fig = plt.figure(figsize= (4,3.5))
axes1 = fig.add_axes([0,0,1,1])


axes1.tick_params(axis='both', which='major', labelsize=16)
##################################################
xmajorLocator   = MultipleLocator(20000)
xmajorFormatter = FormatStrFormatter('%d')
xminorLocator   = MultipleLocator(4000)
axes1.xaxis.set_major_locator(xmajorLocator)
axes1.xaxis.set_major_formatter(xmajorFormatter)
axes1.xaxis.set_minor_locator(xminorLocator)
ymajorLocator   = MultipleLocator(20000)
ymajorFormatter = FormatStrFormatter('%d')
yminorLocator   = MultipleLocator(4000)
axes1.yaxis.set_major_locator(ymajorLocator)
axes1.yaxis.set_major_formatter(ymajorFormatter)
axes1.yaxis.set_minor_locator(yminorLocator)
##################################################
axes1.set_ylabel('Observed Temperature (K)')
axes1.set_xlabel("Model Temperature (K)")
axes1.set_ylim((0, 100000))
axes1.set_xlim((0, 100000))

axes1.plot(x, x, color = 'black', ls = '--')
ax1 = axes1.scatter(desi_T[wd_type == 1], T1[wd_type == 1], s=5, c= StN[wd_type == 1], vmin=0, vmax=30, cmap = cm, zorder = 100)
plt.colorbar(ax1)
axes1.scatter(desi_T[wd_type == 1], T2[wd_type == 1], color = 'red', s=10, alpha = 0.1)
plt.savefig(path + 'fc_T_fit_vs_model_' + model + '.pdf', bbox_inches = 'tight')
plt.show()
plt.close()



fig = plt.figure(figsize= (7,4))
axes1 = fig.add_axes([0,0,1,1])


axes1.tick_params(axis='both', which='major', labelsize=16)
##################################################
xmajorLocator   = MultipleLocator(20000)
xmajorFormatter = FormatStrFormatter('%d')
xminorLocator   = MultipleLocator(4000)
axes1.xaxis.set_major_locator(xmajorLocator)
axes1.xaxis.set_major_formatter(xmajorFormatter)
axes1.xaxis.set_minor_locator(xminorLocator)
#ymajorLocator   = MultipleLocator(0.1)
#ymajorFormatter = FormatStrFormatter('%.1f')
#yminorLocator   = MultipleLocator(0.02)
#axes1.yaxis.set_major_locator(ymajorLocator)
#axes1.yaxis.set_major_formatter(ymajorFormatter)
#axes1.yaxis.set_minor_locator(yminorLocator)
##################################################
axes1.set_ylabel(r'$\frac{T_{obs} - T_{model}}{T_{err}}$')
axes1.set_xlabel("Model Temperature (K)")
#axes1.set_ylim((-0.5, 0.5))
axes1.set_xlim((0, 100000))

ax1 = plt.scatter(desi_T[wd_type == 1], (T1[wd_type == 1] - desi_T[wd_type == 1])/(T1_err[wd_type == 1]), c= StN[wd_type == 1], vmin=0, vmax=30, cmap = cm, s=5, zorder = 100)
plt.colorbar(ax1)
plt.scatter(desi_T[wd_type == 1], (T2[wd_type == 1] - desi_T[wd_type == 1])/(T2_err[wd_type == 1]), color = 'red', s=10, alpha = 0.1)
plt.axhline(0, color = '0.5', ls = '--')

plt.savefig(path + 'fc_errT_fit_vs_model_' + model + '.pdf', bbox_inches = 'tight')
#axes1.set_ylim((-0.25, 0.25))
axes1.set_xlim((0, 100000))
plt.savefig(path + 'fc_errT_fit_vs_model_' + model + '_zoom.pdf', bbox_inches = 'tight')

plt.show()
plt.close()



fig = plt.figure(figsize= (4,4))
axes1 = fig.add_axes([0,0,1,1])
axes1.tick_params(axis='both', which='major', labelsize=16)
##################################################
xmajorLocator   = MultipleLocator(0.5)
xmajorFormatter = FormatStrFormatter('%.1f')
xminorLocator   = MultipleLocator(0.1)
axes1.xaxis.set_major_locator(xmajorLocator)
axes1.xaxis.set_major_formatter(xmajorFormatter)
axes1.xaxis.set_minor_locator(xminorLocator)
ymajorLocator   = MultipleLocator(0.5)
ymajorFormatter = FormatStrFormatter('%.1f')
yminorLocator   = MultipleLocator(0.1)
axes1.yaxis.set_major_locator(ymajorLocator)
axes1.yaxis.set_major_formatter(ymajorFormatter)
axes1.yaxis.set_minor_locator(yminorLocator)
##################################################
axes1.set_ylabel('Observed logg')
axes1.set_xlabel("Model logg")
axes1.set_ylim((3, 10))
axes1.set_xlim((6, 10))
plt.scatter(desi_g[wd_type == 1], g1[wd_type == 1], s=5, c= StN[wd_type == 1], vmin=0, vmax=30, cmap = cm, zorder = 100)
plt.colorbar(ax1)
plt.scatter(desi_g[wd_type == 1], g2[wd_type == 1], color = 'red', s=10, alpha = 0.1)
x = np.arange(6, 11, 1)
plt.plot(x, x, color = 'black', ls = '--')
plt.savefig(path + 'g_fit_vs_model_' + model + '.pdf', bbox_inches = 'tight')
plt.show()
plt.close()


fig = plt.figure(figsize= (7,4))
axes1 = fig.add_axes([0,0,1,1])
axes1.tick_params(axis='both', which='major', labelsize=16)
##################################################
xmajorLocator   = MultipleLocator(0.5)
xmajorFormatter = FormatStrFormatter('%.1f')
xminorLocator   = MultipleLocator(0.1)
axes1.xaxis.set_major_locator(xmajorLocator)
axes1.xaxis.set_major_formatter(xmajorFormatter)
axes1.xaxis.set_minor_locator(xminorLocator)
#ymajorLocator   = MultipleLocator(0.1)
#ymajorFormatter = FormatStrFormatter('%.1f')
#yminorLocator   = MultipleLocator(0.02)
#axes1.yaxis.set_major_locator(ymajorLocator)
#axes1.yaxis.set_major_formatter(ymajorFormatter)
#axes1.yaxis.set_minor_locator(yminorLocator)
##################################################
axes1.set_ylabel(r'$\frac{g_{obs} - g_{model}}{g_{err}}$')
axes1.set_xlabel("Model logg")
axes1.set_xlim((6, 10))
#axes1.set_ylim((-0.5, 0.5))
axes1.scatter(desi_g[wd_type == 1], (g1[wd_type == 1] - desi_g[wd_type == 1])/(g1_err[wd_type == 1]), s=5, c= StN[wd_type == 1], vmin=0, vmax=30, cmap = cm, zorder = 100)
plt.colorbar(ax1)
axes1.scatter(desi_g[wd_type == 1], (g2[wd_type == 1] - desi_g[wd_type == 1])/(g2_err[wd_type == 1]), color = 'red', s=10, alpha = 0.1)
axes1.axhline(0, color = '0.5', ls = '--')
plt.savefig(path + 'errg_fit_vs_model_' + model + '.pdf', bbox_inches = 'tight')
axes1.set_xlim((6, 10))
#axes1.set_ylim((-0.25, 0.25))
plt.savefig(path + 'errg_fit_vs_model_' + model + '_zoom.pdf', bbox_inches = 'tight')
plt.show()
plt.close()




plt.scatter(desi_T, desi_g, s = 1, color = 'black')
plt.scatter(T1, g1, c= StN, vmin=0, vmax=30, s = 1)
plt.colorbar(ax1)
plt.show()
plt.close()