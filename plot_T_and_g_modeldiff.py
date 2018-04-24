import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
model = 'da2014'
path = '/Users/christophermanser/Storage/PhD_files/DESI/2percent/'
f = path + 'WD_list_' + model + '.lst'
x = np.arange(5000, 100000, 1000)
cm = plt.cm.get_cmap('viridis')
name, T1, T1_err, g1, g1_err, T2, T2_err, g2, g2_err, rv, StN = np.genfromtxt(f, unpack = True)
desi_name, desi_T, desi_g, wd_type, mag_desi = np.genfromtxt(path + 'desi_WD_info.dat', unpack = True)


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
plt.savefig(path + 'T_fit_vs_model_' + model + '.pdf', bbox_inches = 'tight')
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
ymajorLocator   = MultipleLocator(0.1)
ymajorFormatter = FormatStrFormatter('%.1f')
yminorLocator   = MultipleLocator(0.02)
axes1.yaxis.set_major_locator(ymajorLocator)
axes1.yaxis.set_major_formatter(ymajorFormatter)
axes1.yaxis.set_minor_locator(yminorLocator)
##################################################
axes1.set_ylabel(r'1 - $\frac{T_{obs}}{T_{model}}$')
axes1.set_xlabel("Model Temperature (K)")
axes1.set_ylim((-0.5, 0.5))
axes1.set_xlim((0, 100000))

ax1 = plt.scatter(desi_T[wd_type == 1], 1 - (T1[wd_type == 1])/(desi_T[wd_type == 1]), c= StN[wd_type == 1], vmin=0, vmax=30, cmap = cm, s=5, zorder = 100)
plt.colorbar(ax1)
plt.scatter(desi_T[wd_type == 1], 1 - (T2[wd_type == 1])/(desi_T[wd_type == 1]), color = 'red', s=10, alpha = 0.1)
plt.axhline(0, color = '0.5', ls = '--')

plt.savefig(path + 'precentT_fit_vs_model_' + model + '.pdf', bbox_inches = 'tight')
axes1.set_ylim((-0.25, 0.25))
axes1.set_xlim((0, 100000))
plt.savefig(path + 'precentT_fit_vs_model_' + model + '_zoom.pdf', bbox_inches = 'tight')

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
ymajorLocator   = MultipleLocator(0.1)
ymajorFormatter = FormatStrFormatter('%.1f')
yminorLocator   = MultipleLocator(0.02)
axes1.yaxis.set_major_locator(ymajorLocator)
axes1.yaxis.set_major_formatter(ymajorFormatter)
axes1.yaxis.set_minor_locator(yminorLocator)
##################################################
axes1.set_ylabel(r'1 - $\frac{g_{obs}}{g_{model}}$')
axes1.set_xlabel("Model logg")
axes1.set_xlim((6, 10))
axes1.set_ylim((-0.5, 0.5))
axes1.scatter(desi_g[wd_type == 1], 1 - (g1[wd_type == 1])/(desi_g[wd_type == 1]), s=5, c= StN[wd_type == 1], vmin=0, vmax=30, cmap = cm, zorder = 100)
plt.colorbar(ax1)
axes1.scatter(desi_g[wd_type == 1], 1 - (g2[wd_type == 1])/(desi_g[wd_type == 1]), color = 'red', s=10, alpha = 0.1)
axes1.axhline(0, color = '0.5', ls = '--')
plt.savefig(path + 'precentg_fit_vs_model_' + model + '.pdf', bbox_inches = 'tight')
axes1.set_xlim((6, 10))
axes1.set_ylim((-0.25, 0.25))
plt.savefig(path + 'precentg_fit_vs_model_' + model + '_zoom.pdf', bbox_inches = 'tight')
plt.show()
plt.close()




plt.scatter(desi_T, desi_g, s = 1, color = 'black')
plt.scatter(T1, g1, c= StN, vmin=0, vmax=30, s = 1)
plt.colorbar(ax1)
plt.show()
plt.close()