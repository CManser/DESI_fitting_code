import numpy as np
import sys
import matplotlib.pyplot as plt
from scipy import optimize
import fitting_scripts
import fitting_scripts
model_c='da2014'
basedir='/Users/christophermanser/Storage/PhD_files/DESI/WDFitting'
c = 299792.458 # Speed of light in km/s
plot = True
# Loads the input spectrum as sys.argv[1], first input, and normalises
spectra=np.loadtxt(sys.argv[1],usecols=(0,1,2),unpack=True).transpose()
spectra = spectra[np.isnan(spectra[:,1])==False & (spectra[:,0]>3500)]
spec_w = spectra[:,0]
spec_n, cont_flux = fitting_scripts.norm_spectra(spectra)
#load lines to fit and crops them
line_crop = np.loadtxt(basedir+'/line_crop.dat')
l_crop = line_crop[(line_crop[:,0]>spec_w.min()) & (line_crop[:,1]<spec_w.max())]
#fit entire grid to find good starting point
best=fitting_scripts.fit_line(spec_n, l_crop,model=model_c)
first_T, first_g = best[2][0], best[2][1]
all_chi, all_TL  = best[4], best[3]

#------find starting point for secondary solution
if first_T < 13000.:
    other_TL, other_chi = all_TL[all_TL[:,0]>13000.], all_chi[all_TL[:,0]>13000.]
    other_sol= other_TL[other_chi==np.min(other_chi)]
elif first_T > 13000.:
    other_TL, other_chi = all_TL[all_TL[:,0]<13000.], all_chi[all_TL[:,0]<13000.]
    other_sol= other_TL[other_chi==np.min(other_chi)]

# Find best fitting model and then calculate error
new_best= optimize.fmin(fitting_scripts.fit_func,(first_T,first_g,10.),
                        args=(spec_n,l_crop,model_c,0),disp=0,xtol=1.,
                        ftol=1.,full_output=1)
other_T = optimize.fmin(fitting_scripts.err_func,(first_T,first_g),
                        args=(new_best[0][2],new_best[1],spec_n,l_crop,model_c),
                        disp=0,xtol=1.,ftol=1.,full_output=0)
best_T, best_g, best_rv = new_best[0][0], new_best[0][1], new_best[0][2]
print("First solution")
print("T = ", best_T, abs(best_T-other_T[0]))
print("logg = ", best_g/100, abs(best_g-other_T[1])/100)
print("rv =",best_rv)
#repeat fit for secondary solution
sec_best = optimize.fmin(fitting_scripts.fit_func,(other_sol[0][0],other_sol[0][1],
                         best_rv),args=(spec_n,l_crop,model_c,0),disp=0,xtol=1.,
                         ftol=1.,full_output=1)
other_T2 = optimize.fmin(fitting_scripts.err_func,(other_sol[0][0],other_sol[0][1]),
                         args=(best_rv,sec_best[1],spec_n,l_crop,model_c),disp=0,
                         xtol=1.,ftol=1.,full_output=0)
s_best_T, s_best_g, s_best_rv = sec_best[0][0], sec_best[0][1], sec_best[0][2]
print("\nSecond solution")
print("T = ", s_best_T, abs(s_best_T-other_T2[0]))
print("logg = ", s_best_g/100, abs(s_best_g-other_T2[1])/100)

#=======================plotting===============================================
if plot == True:
    # Get and save the 2 best lines from the spec and model, and the full models
    lines_s,lines_m,mod_n=fitting_scripts.fit_func((best_T,best_g,best_rv),spec_n,
                                                   l_crop,models=model_c,mode=1)
    lines_s_o,lines_m_o,mod_n_o=fitting_scripts.fit_func((s_best_T,s_best_g,s_best_rv),
                                                         spec_n,l_crop,models=model_c,
                                                         mode=1)
    fig=plt.figure(figsize=(8,5))
    ax1 = plt.subplot2grid((1,4), (0, 3))
    step = 0
    for i in range(1,6): # plots Halpha (i=0) to H6 (i=5)
        min_p   = lines_s[i][:,0][lines_s[i][:,1]==np.min(lines_s[i][:,1])][0]
        min_p_o = lines_s_o[i][:,0][lines_s_o[i][:,1]==np.min(lines_s_o[i][:,1])][0]
        ax1.plot(lines_s[i][:,0]-min_p,lines_s[i][:,1]+step,color='k')
        ax1.plot(lines_s[i][:,0]-min_p,lines_m[i]+step,color='r')
        ax1.plot(lines_s_o[i][:,0]-min_p_o,lines_m_o[i]+step,color='g')
        step+=0.5
    xticks = ax1.xaxis.get_major_ticks()
    ax1.set_xticklabels([])
    ax1.set_yticklabels([])

    ax2 = plt.subplot2grid((1,4), (0, 0),colspan=3)
    ax2.plot(spec_w,spectra[:,1],color='k')
    # Adjust the flux of models to match the spectrum
    mod_n[np.isnan(mod_n)], mod_n_o[np.isnan(mod_n_o)] = 0.0, 0.0
    check_f_spec=spectra[:,1][(spec_w>4500.) & (spec_w<4700.)]
    check_f_model=mod_n[:,1][(mod_n[:,0]>4500.) & (mod_n[:,0]<4700.)]
    adjust=np.average(check_f_model)/np.average(check_f_spec)
    ax2.plot(mod_n[:,0]*(best_rv+c)/c,mod_n[:,1]/adjust,color='r')
    check_f_model_o=mod_n_o[:,1][(mod_n_o[:,0]>4500.) & (mod_n_o[:,0]<4700.)]
    adjust_o=np.average(check_f_model_o)/np.average(check_f_spec)
    ax2.plot(mod_n_o[:,0]*(best_rv+c)/c,mod_n_o[:,1]/adjust_o,color='g')

    ax2.set_ylabel(r'F$_{\lambda}$ [erg cm$^{-2}$ s$^{-1} \AA^{-1}$]',fontsize=12)
    ax2.set_xlabel(r'Wavelength $(\AA)$',fontsize=12)
    plt.show()
    plt.close()
