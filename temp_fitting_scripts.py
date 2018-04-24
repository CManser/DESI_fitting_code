import numpy as np
basedir='/Users/christophermanser/Storage/PhD_files/DESI/WDFitting'
c = 299792.458 # Speed of light in km/s

def air_to_vac(wavin):
    """Converts spectra from air wavelength to vacuum to compare to models"""
    return wavin/(1.0 + 2.735182e-4 + 131.4182/wavin**2 + 2.76249e8/wavin**4)


def fit_DA(spectra,inp = 'None', model_c = 'da2014', plot = False):
    import matplotlib.pyplot as plt
    from scipy import optimize
    # Loads the input spectrum as inp (read in spectrum), first input, and normalises
    #load lines to fit and crops them
    spec_w = spectra[:,0]
    #spectra[:,2]=spectra[:,2]**-0.5
    spec_n, cont_flux = norm_spectra(spectra)
    line_crop = np.loadtxt(basedir+'/line_crop.dat')
    l_crop = line_crop[(line_crop[:,0]>spec_w.min()) & (line_crop[:,1]<spec_w.max())]
    #fit entire grid to find good starting point
    best=fit_line(spec_n, l_crop,model=model_c)
    first_T, first_g = best[2][0], best[2][1]
    all_chi, all_TL  = best[4], best[3]

    #------find starting point for secondary solution
    if first_T <= 13000.:
        other_TL, other_chi = all_TL[all_TL[:,0]>=13000.], all_chi[all_TL[:,0]>=13000.]
        other_sol = other_TL[other_chi==np.min(other_chi)]
    elif first_T > 13000.:
        other_TL, other_chi = all_TL[all_TL[:,0]<13000.], all_chi[all_TL[:,0]<13000.]
        other_sol = other_TL[other_chi==np.min(other_chi)]

    # Find best fitting model and then calculate error
    new_best= optimize.fmin(fit_func,(first_T,first_g,10.),args=(spec_n,l_crop,
                            model_c,0),disp=0,xtol=1.,ftol=1.,full_output=1)
    other_T = optimize.fmin(err_func,(first_T,first_g),args=(new_best[0][2],
                            new_best[1],spec_n,l_crop,model_c),disp=0,xtol=1.,
                            ftol=1.,full_output=0)
    best_T, best_g, best_rv = new_best[0][0], new_best[0][1], new_best[0][2]
    print("\nFirst solution")
    print("T = ", best_T, abs(best_T-other_T[0]))
    print("logg = ", best_g/100, abs(best_g-other_T[1])/100)
    print("rv =",best_rv)
    #repeat fit for secondary solution
    sec_best = optimize.fmin(fit_func,(other_sol[0][0],other_sol[0][1],best_rv),
                             args=(spec_n,l_crop,model_c,0),disp=0,xtol=1.,ftol=1.,
                             full_output=1)
    other_T2 = optimize.fmin(err_func,(other_sol[0][0],other_sol[0][1]),
                             args=(best_rv,sec_best[1],spec_n,l_crop,model_c),disp=0,
                             xtol=1.,ftol=1.,full_output=0)
    s_best_T, s_best_g, s_best_rv = sec_best[0][0], sec_best[0][1], sec_best[0][2]
    print("\nSecond solution")
    print("T = ", s_best_T, abs(s_best_T-other_T2[0]))
    print("logg = ", s_best_g/100, abs(s_best_g-other_T2[1])/100)

    StN = np.sum( spectra[:,1][(spec_w >= 4500.0) & (spec_w <= 4750.0)] / spectra[:,2][(spec_w >= 4500.0) & (spec_w <= 4750.0)] ) / spectra[:,1][(spec_w >= 4500.0) & (spec_w <= 4750.0)].size
    print('StN = %.1f'%(StN))
    #=======================plotting===============================================
    if plot == True:
        # Get and save the 2 best lines from the spec and model, and the full models
        lines_s,lines_m,mod_n=fit_func((best_T,best_g,best_rv),spec_n,l_crop,
                                       models=model_c,mode=1)
        lines_s_o,lines_m_o,mod_n_o=fit_func((s_best_T,s_best_g,s_best_rv),spec_n,
                                             l_crop,models=model_c,mode=1)
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
        ax2.plot(mod_n_o[:,0]*(best_rv+c)/c,mod_n_o[:,1]/adjust_o,color='g', ls = '--')

        ax2.set_ylabel(r'F$_{\lambda}$',fontsize=12)
        ax2.set_xlabel(r'Wavelength $(\AA)$',fontsize=12)
    
        ax2.set_xlim(3500, 6000)
        ax2.set_ylim([spectra[:,1].min(), spectra[:,1].max()])
    
        ax2.text(4500, spectra[:,1].min() + 0.96*( spectra[:,1].max() - spectra[:,1].min()), 'StN = %.1f'%(StN), bbox={'edgecolor':'white', 'facecolor':'white', 'alpha':1, 'pad':2}, fontsize = 8)
        name, T_desi, g_desi, wd_type, mag_desi = np.genfromtxt('/Users/christophermanser/Storage/PhD_files/DESI/2percent/desi_wd_info.dat', unpack = True)
        tmp = int(inp.split('/')[-1].split('_')[0])
        mag = mag_desi[np.where(tmp == name)]
        ax2.text(  5200, spectra[:,1].min() + 0.96*( spectra[:,1].max() - spectra[:,1].min()), 'Mag = %.2f'%(mag), bbox={'edgecolor':'white', 'facecolor':'white', 'alpha':1, 'pad':2}, fontsize = 8)
        ax2.text(  4500, spectra[:,1].min() + 0.92*( spectra[:,1].max() - spectra[:,1].min()), 'T = %d +/- %d'%(best_T, abs(best_T-other_T[0])), bbox={'edgecolor':'white', 'facecolor':'white', 'alpha':1, 'pad':2}, fontsize = 8)
        ax2.text(  5250, spectra[:,1].min() + 0.92*( spectra[:,1].max() - spectra[:,1].min()), 'logg = %.2f +/- %.2f'%(best_g/100, abs(best_g-other_T[1])/100), bbox={'edgecolor':'white', 'facecolor':'white', 'alpha':1, 'pad':2}, fontsize = 8)
        ax2.text(  4500, spectra[:,1].min() + 0.88*( spectra[:,1].max() - spectra[:,1].min()), 'T2 = %d +/- %d'%(s_best_T, abs(s_best_T-other_T2[0])), bbox={'edgecolor':'white', 'facecolor':'white', 'alpha':1, 'pad':2}, fontsize = 8)
        ax2.text(  5250, spectra[:,1].min() + 0.88*( spectra[:,1].max() - spectra[:,1].min()), 'logg2 = %.2f +/- %.2f'%(s_best_g/100, abs(s_best_g-other_T2[1])/100), bbox={'edgecolor':'white', 'facecolor':'white', 'alpha':1, 'pad':2}, fontsize = 8) 
        ax2.text(  4500, spectra[:,1].min() + 0.84*( spectra[:,1].max() - spectra[:,1].min()), 'T DESI = %d'%(T_desi[np.where(tmp == name)]), bbox={'edgecolor':'white', 'facecolor':'white', 'alpha':1, 'pad':2}, fontsize = 8)
        ax2.text(  5250, spectra[:,1].min() + 0.84*( spectra[:,1].max() - spectra[:,1].min()), 'logg DESI = %.2f'%(g_desi[np.where(tmp == name)]), bbox={'edgecolor':'white', 'facecolor':'white', 'alpha':1, 'pad':2}, fontsize = 8)    
        plt.savefig('/Users/christophermanser/Storage/PhD_files/DESI/WD_ascii/WD_fits/20171117/' + inp[len(basedir):-4] + '-' + model_c + '.png', dpi = 300, bbox_inches = 'tight')
        plt.close()
    return best_T, best_g, best_rv


def err_func(x,rv,valore,specn,lcrop,models='da2014'):
    """Script finds errors by minimising function at chi+1 rather than chi
       Requires: x; rv - initial guess of T, g; rv
       valore - the chi value of the best fit
       specn/lcrop - normalised spectrum / list of cropped lines to fit"""
    tmp = tmp_func(x[0], x[1], rv, specn, lcrop, models)
    if tmp != 1: return abs(tmp[3]-(valore+1.)) #this is quantity that gets minimized 
    else: return 1E30


def fit_func(x,specn,lcrop,models='da2014',mode=0):
    """Requires: x - initial guess of T, g, and rv
       specn/lcrop - normalised spectrum / list of cropped lines to fit
       mode=0 is for finding bestfit, mode=1 for fitting & retriving specific model """
    tmp = tmp_func(x[0], x[1], x[2], specn, lcrop, models)
    if tmp == 1: return 1E30
    elif mode==0: return tmp[3] #this is the quantity that gets minimized
    elif mode==1: return tmp[0], tmp[1], tmp[2]


def fit_line(_sn, l_crop, model='da2014'):
    """Input _sn, l_crop <- Normalised spectrum & a cropped line list
    Optional:
        model='da2014'  : 'da2014' or 'pier'
    Calc and return chi2, list of arrays of spectra, and scaled models at lines """
    from scipy import interpolate
    #load normalised models and linearly interp models onto spectrum wave
    m_wave,m_flux_n,m_param = norm_models(model=model)
    sn_w = _sn[:,0]
    m_flux_n_i = interpolate.interp1d(m_wave,m_flux_n,kind='linear')(sn_w)
    #Crops models and spectra in a line region, renorms models, calculates chi2
    tmp_lines_m, lines_s, l_chi2 = [],[],[]
    for i in range(len(l_crop)):
        l_c0,l_c1 = l_crop[i,0],l_crop[i,1] 
        l_m = m_flux_n_i.transpose()[(sn_w>=l_c0)&(sn_w<=l_c1)].transpose()
        l_s = _sn[(sn_w>=l_c0)&(sn_w<=l_c1)]
        l_m = l_m*np.sum(l_s[:,1])/np.sum(l_m,axis=1).reshape([len(l_m),1])
        l_chi2.append( np.sum(((l_s[:,1]-l_m)/l_s[:,2])**2,axis=1) )
        tmp_lines_m.append(l_m)
        lines_s.append(l_s)
    #mean chi2 over lines and stores best model lines for output
    lines_chi2, lines_m = np.sum(np.array(l_chi2),axis=0), []
    is_best = lines_chi2==lines_chi2.min()
    for i in range(len(l_crop)): lines_m.append(tmp_lines_m[i][is_best][0])
    best_TL = m_param[is_best][0]
    return  lines_s,lines_m,best_TL,m_param,lines_chi2


def interpolating_model_DA(temp,grav,m_type='da2014'):
    """Interpolate model atmospheres given an input Teff and logg
    models are saved as a numpy array to increase speed"""
    # PARAMETERS # 
    dir_models = basedir + '/WDModels_Koester.'+m_type+'_npy/'
    if m_type=="pier":
        teff=np.array([1500.,1750.,2000.,2250.,2500.,2750.,3000.,3250.,3500.,
                       3750.,4000.,4250.,4500.,4750.,5000.,5250.,5500.,6000.,
                       6500.,7000.,7500.,8000.,8500.,9000.,9500.,10000.,10500.,
                       11000.,11500.,12000.,12500.,13000.,13500.,14000.,14500.,
                       15000.,15500.,16000.,16500.,17000.,20000.,25000.,30000.,
                       35000.,40000.,45000.,50000.,55000.,60000.,65000.,70000.,
                       75000.,80000.,85000.,90000.])
        logg=np.array([6.50,7.00,7.50,7.75,8.00,8.25,8.50,9.00,9.50])
    elif m_type=="da2014":
        teff=np.array([6000.,6250.,6500.,6750.,7000.,7250.,7500.,7750.,8000.,
                       8250.,8500.,8750.,9000.,9250.,9500.,9750.,10000.,10100.,
                       10200.,10250.,10300.,10400.,10500.,10600.,10700.,10750.,
                       10800.,10900.,11000.,11100.,11200.,11250.,11300.,11400.,
                       11500.,11600.,11700.,11750.,11800.,11900.,12000.,12100.,
                       12200.,12250.,12300.,12400.,12500.,12600.,12700.,12750.,
                       12800.,12900.,13000.,13500.,14000.,14250.,14500.,14750.,
                       15000.,15250.,15500.,15750.,16000.,16250.,16500.,16750.,
                       17000.,17250.,17500.,17750.,18000.,18250.,18500.,18750.,
                       19000.,19250.,19500.,19750.,20000.,21000.,22000.,23000.,
                       24000.,25000.,26000.,27000.,28000.,29000.,30000.,35000.,
                       40000.,45000.,50000.,55000.,60000.,65000.,70000.,75000.,
                       80000.,90000.,100000.])
        logg=np.array([4.00,4.25,4.50,4.75,5.00,5.25,5.50,5.75,6.00,6.25,6.50,
                       6.75,7.00,7.25,7.50,7.75,8.00,8.25,8.50,8.75,9.00,9.25,
                       9.50])
    if (m_type=='pier') & (temp<1500. or temp>90000. or grav<6.5 or grav>9.5): 
        return [],[]
    elif (m_type=='da2014') & (temp<6000. or temp>100000. or grav<4.0 or grav>9.5): 
        return [],[]
    # INTERPOLATION #
    g1,g2 = np.max(logg[logg<=grav]),np.min(logg[logg>=grav])
    if g1!=g2: g = (grav-g1)/(g2-g1)
    else: g=0
    t1,t2 = np.max(teff[teff<=temp]),np.min(teff[teff>=temp])
    if t1!=t2: t = (temp-t1)/(t2-t1)          
    else: t=0	
    if m_type =='da2014': models = ['da%06d_%d_2.7.npy'%(i, j*100) for i in [t1,t2] 
                                    for j in [g1,g2]]
    else: models = ['WD_%.2f_%d.0.npy'%(j, i) for i in [t1,t2] for j in [g1,g2]]
    try:
        m11, m12 = np.load(dir_models+models[0]), np.load(dir_models+models[1])	
        m21, m22 = np.load(dir_models+models[2]), np.load(dir_models+models[3])	
        flux_i = (1-t)*(1-g)*m11[:,1]+t*(1-g)*m21[:,1]+t*g*m22[:,1]+(1-t)*g*m12[:,1]
        return np.dstack((m11[:,0], flux_i))[0]
    except: return [],[]


def norm_models(model='da2014'):
    """ Import Normalised WD Models
    Optional arguments:
        model='da2014': Which model grid to use: List shown below in mode_list
    Return [out_m_wave,norm_m_flux,model_param] """
    model_list = ['da2014','pier']
    if model not in model_list: raise wdfitError('Unknown "model" in norm_models')
    fn, d = '/wdfit.'+model+'.lst', '/WDModels_Koester.'+model+'_npy/'
    model_list  = np.loadtxt(basedir+fn, usecols=[0], dtype=np.string_).astype(str)
    model_param = np.loadtxt(basedir+fn, usecols=[1,2])
    m_spec = np.load(basedir+d+model_list[0])
    m_wave = m_spec[:,0]
    out_m_wave = m_wave[(m_wave>=3400)&(m_wave<=13000)]
    norm_m_flux = np.load(basedir+'/norm_m_flux.'+model+'.npy')
    #
    if out_m_wave.shape[0] != norm_m_flux.shape[1]:
        raise wdfitError('l and f arrays not correct shape in norm_models')
    return [out_m_wave,norm_m_flux,model_param]


def norm_spectra(spectra, add_infinity=True):
    """ Normalised spectra by DA WD continuum regions. Spec of form 
    array([wav,flux,err]) (err not necessary) only works on SDSS spectra region
    Optional:
        Edit n_range_s to change whether region[j] is fitted for a peak or mean'd
        add_infinity=False : add a spline point at [inf,0]
    returns spectra, cont_flux """
    from scipy import interpolate
    start_n=np.array([3770.,3796.,3835.,3895.,3995.,4130.,4490.,4620.,5070.,5200.,
                      6000.,7000.,7550.,8400.])
    end_n=np.array([3795.,3830.,3885.,3960.,4075.,4290.,4570.,4670.,5100.,5300.,
                    6100.,7050.,7600.,8450.])
    n_range_s=np.array(['P','P','P','P','P','P','M','M','M','M','M','M','M','M'])
    if len(spectra[0])>2:
        snr = np.zeros([len(start_n),3])
        spectra[:,2][spectra[:,2]==0.] = spectra[:,2].max()
    else: 
        snr = np.zeros([len(start_n),2])
    wav = spectra[:,0]
    for j in range(len(start_n)):
        if (start_n[j] < wav.max()) & (end_n[j] > wav.min()):
            _s = spectra[(wav>=start_n[j])&(wav<=end_n[j])]
            _w = _s[:,0]
            #Avoids gappy spectra
            k=3 # Check if there are more points than 3
            if len(_s)>k:
                #interpolate onto 10* resolution
                l = np.linspace(_w.min(),_w.max(),(len(_s)-1)*10+1)
                if len(spectra[0])>2:
                    tck = interpolate.splrep(_w,_s[:,1],w=1/_s[:,2], s=2000)
                    #median errors for max/mid point
                    snr[j,2] = np.median(_s[:,2]) / np.sqrt(len(_w))
                else: tck = interpolate.splrep(_w,_s[:,1],s=0.0)
                f = interpolate.splev(l,tck)
                #find maxima and save
                if n_range_s[j]=='P': snr[j,0], snr[j,1] = l[f==f.max()][0], f.max()
                #find mean and save
                elif n_range_s[j]=='M': snr[j,0:2] = np.mean(l), np.mean(f)
                else: print('Unknown n_range_s, ignoring')
    snr = snr[ snr[:,0] != 0 ]
    #t parameter chosen by eye. Position of knots.
    if snr[:,0].max() < 6460: knots = [3000,4900,4100,4340,4860,int(snr[:,0].max()-5)]
    else: knots = [3885,4340,4900,6460]
    if snr[:,0].min() > 3885:
        print('Warning: knots used for spline norm unsuitable for high order fitting')
        knots=knots[1:]
    if (snr[:,0].min() > 4340) or (snr[:,0].max() < 4901): 
        knots=None # 'Warning: knots used probably bad'
    if add_infinity: # Adds points at inf & 0 for spline to fit to err = mean(spec err)
        if snr.shape[1] > 2:
            mean_snr = np.mean(snr[:,2])
            snr = np.vstack([ snr, np.array([90000. ,0., mean_snr ]) ])
            snr = np.vstack([ snr, np.array([100000.,0., mean_snr ]) ])
        else:
            snr = np.vstack([ snr, np.array([90000.,0.]) ])
            snr = np.vstack([ snr, np.array([100000.,0.]) ])
    try: #weight by errors
        if len(spectra[0])>2: 
            tck = interpolate.splrep(snr[:,0],snr[:,1], w=1/snr[:,2], t=knots, k=3)
        else: tck = interpolate.splrep(snr[:,0],snr[:,1], t=knots, k=3)
    except ValueError:
        knots=None
        if len(spectra[0])>2: 
            tck = interpolate.splrep(snr[:,0],snr[:,1], w=1/snr[:,2], t=knots, k=3)
        else: tck = interpolate.splrep(snr[:,0],snr[:,1], t=knots, k=3)
    cont_flux = interpolate.splev(wav,tck).reshape(wav.size, 1)
    spectra_ret = np.copy(spectra)
    spectra_ret[:,1:] = spectra_ret[:,1:]/cont_flux
    return spectra_ret, cont_flux
    

def tmp_func(_T, _g, _rv, _sn, _l, _m):
    from scipy import interpolate
    model=interpolating_model_DA(_T,(_g/100),m_type=_m)
    try: norm_model, m_cont_flux=norm_spectra(model)
    except:
        print("Could not load the model")
        return 1
    else:
        #interpolate normalised model and spectra onto same wavelength scale
        m_wave_n, m_flux_n, sn_w = norm_model[:,0]*(_rv+c)/c, norm_model[:,1], _sn[:,0]
        m_flux_n_i = interpolate.interp1d(m_wave_n,m_flux_n,kind='linear')(sn_w)
        #Initialise: normalised models and spectra in line region, and chi2
        lines_m, lines_s, sum_l_chi2 = [],[],0
        for i in range(len(_l)):
            # Crop model and spec to line
            l_c0,l_c1 = _l[i,0],_l[i,1]
            l_m = m_flux_n_i.transpose()[(sn_w>=l_c0)&(sn_w<=l_c1)].transpose()
            l_s = _sn[(sn_w>=l_c0)&(sn_w<=l_c1)]
            #renormalise models to spectra in line region & calculate chi2+sum
            l_m = l_m*np.sum(l_s[:,1])/np.sum(l_m)
            sum_l_chi2 += np.sum(((l_s[:,1]-l_m)/l_s[:,2])**2)
            lines_m.append(l_m), lines_s.append(l_s)
        return lines_s, lines_m, model, sum_l_chi2


def flux_calib(spec, inp, knot_num = 40):
    """ Takes in spectrum and path to the spectrum. knot_num determines the num
    of knots used in spline fit of response curve. itera = num of interations.
    
    Returns the response function calculated, as well as the best fitting model
    with the temperature and logg"""
    from scipy.interpolate import splev, splrep, interp1d
    import matplotlib.pyplot as plt
    # BLUE
    temp, logg, rv = fit_DA(spec, inp, model_c = 'da2014', plot = False)
    model = interpolating_model_DA(temp,logg/100,m_type='da2014')
    m1, m2 = model[:,0]*(rv+c)/c, model[:,1]
    best_model = np.stack((m1,m2,np.ones(m1.size)), axis=-1)
    
    
    # Takes the model and divides this out of the spec to determine the resp. curve.
    model_flux = interp1d(m1,m2,kind='linear')(spec[:,0])
    response_curve = spec[:,1]/model_flux
    r_tmp_b = response_curve[:]
    wav_min, wav_max = min(spec[:,0]), max(spec[:,0])
    # Spline fit the response curve to get a smooth function.
    tmp = (wav_max - wav_min)/knot_num
    knots = np.linspace(wav_min+tmp, wav_max-tmp, knot_num)
    tck = splrep(spec[:,0], r_tmp_b, w = 1/spec[:,2], t=knots, k = 3)
    p = splev(spec[:,0], tck)
    p /= np.mean(p)
    r_tmp_b /= np.mean(r_tmp_b)
    # Calculate the accurate flux/count ratio
    print('BLUE ARM: Include flux/count calculation')
    rcf_b = p
    
    #RED
    inp_r = inp[:-20] + 'r' + inp[-19:]
    s_r = np.loadtxt(inp_r,usecols=(0,1,2),unpack=True).transpose()
    s_r = s_r[np.isnan(s_r[:,1])==False & (s_r[:,0]>3500)]
    s_r1, s_r2, s_r3 = s_r[:,0].copy(), s_r[:,1].copy(), s_r[:,2].copy()
    
    model_flux_r = interp1d(m1,m2,kind='linear')(s_r1)
    response_curve = s_r2/model_flux_r
    r_tmp_r = response_curve[:]
    wav_min, wav_max = min(s_r1), max(s_r1)
    # Spline fit the response curve to get a smooth function.
    tmp = (wav_max - wav_min)/knot_num
    knots = np.linspace(wav_min+tmp, wav_max-tmp, knot_num)
    tck = splrep(s_r1, r_tmp_r, w = 1/s_r3, t=knots, k = 3)
    p = splev(s_r1, tck)
    p /= np.mean(p)
    r_tmp_r /= np.mean(r_tmp_r)
    # Calculate the accurate flux/count ratio
    print('RED ARM: Include flux/count calculation')
    rcf_r = p
    
    #ZED
    inp_z = inp[:-20] + 'z' + inp[-19:]
    s_z = np.loadtxt(inp_z,usecols=(0,1,2),unpack=True).transpose()
    s_z = s_z[np.isnan(s_z[:,1])==False & (s_z[:,0]>3500)]
    s_z1, s_z2, s_z3 = s_z[:,0].copy(), s_z[:,1].copy(), s_z[:,2].copy()
    
    model_flux_z = interp1d(m1,m2,kind='linear')(s_z1)
    
    response_curve = s_z2/model_flux_z
    r_tmp_z = response_curve[:]
    wav_min, wav_max = min(s_z1), max(s_z1)
    # Spline fit the response curve to get a smooth function.
    tmp = (wav_max - wav_min)/knot_num
    knots = np.linspace(wav_min+tmp, wav_max-tmp, knot_num)
    tck = splrep(s_z1, r_tmp_z, w = 1/s_z3, t=knots, k = 3)
    p = splev(s_z1, tck)
    p /= np.mean(p)
    r_tmp_z /= np.mean(r_tmp_z)
    # Calculate the accurate flux/count ratio on first itr
    print('ZED ARM: Include flux/count calculation')
    rcf_z = p
    
    
    return [rcf_b, r_tmp_b], [rcf_r, r_tmp_r], [rcf_z, r_tmp_z], [temp, logg, best_model]