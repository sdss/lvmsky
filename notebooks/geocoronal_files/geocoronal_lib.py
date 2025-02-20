# load packages
import numpy as np
from astropy.table import Table
from astropy.time import Time

from astropy.io import fits
from astropy.modeling import models, fitting
from scipy.optimize import least_squares
from collections import Counter

from astropy.coordinates import SkyCoord  # High-level coordinates
import astropy.units as u
import os
import pickle

def load_sky_spectra(t):
    rsshdu = fits.open(t)
    rsshdu.info

    hdr = rsshdu['PRIMARY'].header

    wave=rsshdu['WAVE'].data 

    r1 = rsshdu['FLUX'].data
    r1_ivar = rsshdu['IVAR'].data
    mask = rsshdu['MASK'].data
    sky = rsshdu['SKY'].data
    tab = Table(rsshdu['SLITMAP'].data)

    r1[mask == 1] = np.nan

    skyW = (tab['targettype'] == 'SKY') & (tab['telescope']=='SkyW') & (tab['fibstatus']==0)
    skyE = (tab['targettype'] == 'SKY') & (tab['telescope']=='SkyE') & (tab['fibstatus']==0)
    sci  = (tab['targettype'] == 'SCI') & (tab['fibstatus']==0)

    east_sky_spec = np.nanmedian(r1[skyE,:]+sky[skyE,:],axis=0)
    west_sky_spec = np.nanmedian(r1[skyW,:]+sky[skyW,:],axis=0)

    east_sky_spec_ivar = np.nanmedian(r1_ivar[skyE,:],axis=0)
    west_sky_spec_ivar = np.nanmedian(r1_ivar[skyW,:],axis=0)    

    # fix everything to zenith here!
#    txt = np.genfromtxt(os.getenv("LVMCORE_DIR") + "/etc/lco_extinction.txt")
    txt = np.genfromtxt("lco_extinction.txt")

    lext, ext = txt[:, 0], txt[:, 1]
    ext = np.interp(wave, lext, ext)
    sci_secz = hdr["TESCIAM"]
    skye_secz = hdr["TESKYEAM"]
    skyw_secz = hdr["TESKYWAM"]
    
    east_sky_spec = east_sky_spec / 10 ** (0.4 * ext * (sci_secz)) * 10 ** (0.4 * ext * (skye_secz))
    west_sky_spec = west_sky_spec / 10 ** (0.4 * ext * (sci_secz)) * 10 ** (0.4 * ext * (skyw_secz))
    
    rsshdu.close()

    return wave, west_sky_spec, east_sky_spec, west_sky_spec_ivar, east_sky_spec_ivar, hdr
    
def kin_fit(wave,r1,r1_ivar):
    mm = 6562.79
    wl_ha = (wave > 6559) & (wave < 6566)

    iis_ha = np.where(wl_ha)[0]
    an_amplitude = r1[iis_ha].max()
    an_mean = mm
    an_stddev = 0.65
    g_init = models.Gaussian1D(amplitude=an_amplitude, mean=an_mean, stddev=an_stddev)+models.Const1D(0)

    fit_g = fitting.LevMarLSQFitter()

    g = fit_g(g_init, wave[iis_ha], r1[iis_ha], weights=r1_ivar[iis_ha])

    LVM_FIBER_DIAMETER = 35.3 # arcsec
    fib_area = np.pi*(LVM_FIBER_DIAMETER/2.)**2

    flux = g.amplitude_0.value * np.sqrt(2.*np.pi*g.stddev_0.value**2)/fib_area #erg/s/cm^2/arcsec^2
    mean = g.mean_0.value
    stddev = g.stddev_0.value
 
    try:
        cov_diag = np.diag(fit_g.fit_info['param_cov'])
        flux_err = (g.amplitude_0.value + np.sqrt(cov_diag[0])) * np.sqrt(2.*np.pi*g.stddev_0.value**2)/fib_area - flux
    except:
#        print(f"covariances didn't work, line: {line}, flux: {flux}")
        flux_err = np.nan
        

    return flux, mean, stddev, g(wave), flux_err


def build_geocoronal_database(DIR_redux, ver, rebuild_all=True,verbose=False):
# point this at wherever you want to do the geocoronal fits
# ver ='1.1.1'             # which DRP version?
# DIR_redux = '/data/LVM/' # What is the base directory for the LVM data?
# rebuild_all = True       # start from an existing fits table or rebuild the whole thing?
# verbose = True           # print out some info on progress while going through the drpall files
# geocoronal_lib.build_geocoronal_database('/data/LVM/','1.1.1',True)

    drpall = Table.read(DIR_redux+'sdsswork/lvm/spectro/redux/'+ver+'/drpall-'+ver+'.fits')

    if rebuild_all:
        drpall_geocoronal = drpall
        drpall_geocoronal['geo_flux'] = 0.
        drpall_geocoronal['geo_fluxerr'] = 0.
        drpall_geocoronal['geo_whichsky'] = " " * 20
        drpall_geocoronal['geo_WHAM'] = " " * 20
        drpall_geocoronal['geo_baseline'] = 0.
        drpall_geocoronal['geo_sh_hght'] = 0.

        # shadow height in the drpall file seems to be incorrect, 
        # this code  fixes it but takes a little bit of time to run

        if True:
            from shadow_height_lib import shadow_calc
            from astropy.time import Time

            # initialize the shadow height calculation
            s=shadow_calc()

            for dd in drpall_geocoronal:
                time = Time(dd['obstime'],format='isot', scale='utc')
                jj = time.jd

                cc = SkyCoord(dd['skye_ra'],dd['skye_dec'],frame='icrs',unit=u.deg)
                s.set_coordinates(cc.ra.deg,cc.dec.deg)
                s.update_time(jd=jj)
                dd['skye_sh_hght'] = s.get_heights(return_heights=True, unit="km")

                cc = SkyCoord(dd['skyw_ra'],dd['skyw_dec'],frame='icrs',unit=u.deg)
                s.set_coordinates(cc.ra.deg,cc.dec.deg)
                s.update_time(jd=jj)
                dd['skyw_sh_hght'] = s.get_heights(return_heights=True, unit="km")

                cc = SkyCoord(dd['sci_ra'],dd['sci_dec'],frame='icrs',unit=u.deg)
                s.set_coordinates(cc.ra.deg,cc.dec.deg)
                s.update_time(jd=jj)
                dd['sci_sh_hght'] = s.get_heights(return_heights=True, unit="km")

    else:
        drpall_geocoronal = Table.read('drpall-'+ver+'_geo.fits')

    new = (drpall_geocoronal['exptime']==900) & (drpall_geocoronal['geo_flux']==0) # & (drpall_geocoronal['moon_alt']<0) 
    geocoronal_data = [DIR_redux+dd['location'] for dd in drpall_geocoronal[new]]

    # load in the WHAM darkest field list
    wham = Table.read('darkest_WHAM_list.csv',format='csv')
    cc_wham = SkyCoord(wham['ra[deg]'],wham['dec[deg]'],frame='icrs',unit=u.deg)
    
    
    for t in geocoronal_data:
        if verbose: # print some stuff
            if (np.sum(drpall_geocoronal['geo_flux']!=0) % 100 == 0):
                progress = np.sum(drpall_geocoronal['geo_flux']!=0)
                total = len(geocoronal_data)
                print(f'{progress} out of {total}')
        ii=t.find('sdsswork')

        dd = drpall_geocoronal[t[ii:] == drpall_geocoronal['location']]

        cc_e = SkyCoord(dd['skye_ra'],dd['skye_dec'],frame='icrs',unit=u.deg)
        cc_w = SkyCoord(dd['skyw_ra'],dd['skyw_dec'],frame='icrs',unit=u.deg)

        # find the WHAM field
        if np.sum(cc_e.separation(cc_wham).deg < 1):
            cc = cc_e
            which_sky = 'SKYE'
            wham_name = wham['name'][cc.separation(cc_wham).deg < 1][0]
        elif np.sum(cc_w.separation(cc_wham).deg < 1):
            cc = cc_w
            which_sky = 'SKYW'
            wham_name = wham['name'][cc.separation(cc_wham).deg < 1][0]
        else: # don't bother with fields without WHAM sky
            continue

        if os.path.exists(t):
            wave, swm, sem, swm_ivar, sem_ivar, hdr = load_sky_spectra(t)
        else:
            continue

        if (which_sky == 'SKYE'):
            # find the WHAM field
            sm = sem
            sm_ivar = sem_ivar
        elif (which_sky == 'SKYW'):
            sm = swm
            sm_ivar = swm_ivar


        sub=(wave > 6530) & (wave < 6600)
        # test if baseline is bad
        baseline = np.median(sm[sub])
        if (baseline < 0) | (baseline > 4e-14) | (np.median(sm[(wave > 6530) & (wave < 6540)]) < 0):
            continue

        iis_ha = (wave > 6561) & (wave < 6565)


        if (np.sum(sm[iis_ha]) > np.median(sm[sub])*np.sum(iis_ha)) & (np.median(sm[sub]) > 0) :

            # fit the spectrum
            flux, mean, stddev, fit, flux_err = kin_fit(wave[sub],sm[sub],sm_ivar[sub])

            drpall_geocoronal['geo_flux'][t[ii:] == drpall_geocoronal['location']] = flux
            drpall_geocoronal['geo_fluxerr'][t[ii:] == drpall_geocoronal['location']] = flux_err
            drpall_geocoronal['geo_baseline'][t[ii:] == drpall_geocoronal['location']] = baseline
            drpall_geocoronal['geo_whichsky'][t[ii:] == drpall_geocoronal['location']] = which_sky
            drpall_geocoronal['geo_WHAM'][t[ii:] == drpall_geocoronal['location']] = wham_name
            drpall_geocoronal['geo_sh_hght'][t[ii:] == drpall_geocoronal['location']] = dd[which_sky.lower()+'_sh_hght']

    drpall_geocoronal.write('drpall-'+ver+'_geo.fits',overwrite=True)

# Define the models

# tau is fixed across all data, r0 and rbgr are fixed per night
def model(params, xx, nights):
    x = xx/1e3
    tau = params[0]  # First parameter is tau
    r0 = params[1:len(np.unique(nights)) + 1]  # Normalization for each night
    rbgr = params[len(np.unique(nights)) + 1:]  # Background for each night
    # Calculate the exponential function for each night
    return r0[nights.astype(int)] * np.exp(-x / tau) + rbgr[nights.astype(int)]
# Define the residual function
def residuals_errs(params, x, y, nights,errors):
    return (model(params, x, nights) - y) / errors

def model_fixrbgr(params, xx, tau,rbgr):
    # fix rbgr
    x = xx/1e3
    r0 = params  # Normalization for each night
    # Calculate the exponential function for each night
    return r0 * np.exp(-x / tau) + rbgr
# Define the residual function
def residuals_errs_fixrbgr(params, x, y, tau,rbgr,errors):
    return (model_fixrbgr(params, x, tau,rbgr) - y) / errors

def model_global(params, xx):
    # global fit, no dependence on night
    x = xx/1e3
    tau = params[0]  # First parameter is tau
    r0 = params[1]  # Normalization for each night
    rbgr = params[2]  # Background for each night
    # Calculate the exponential function for each night
    return r0 * np.exp(-x / tau) + rbgr    
# Define the residual function
def residuals_global_errs(params, x, y,errors):
    return (model_global(params, x) - y) / errors

# convert to and from Rayleighs
def cgs_to_R(x):
    return x / 5.661e-18 # returns Rayleigh

def R_to_cgs(x):
    return x * 5.661e-18 # returns erg/s/cm^2/arcsec^2
    
def create_geocoronal_global_model(ver):
# for a given DRP version, return the parameters for a geocoronal model, taking all data together
# returns the parameters of the fit
    drpall_geo = Table.read('drpall-'+ver+'_geo.fits')

    obsT = Time(drpall_geo['obstime'],format='isot', scale='utc') 
    time = obsT.jd % 1 * 24 
    time[time>0]-=24

    iis0 =  (drpall_geo['moon_alt'] < 0) & (drpall_geo['geo_baseline'] < 7e-15) & (drpall_geo['geo_baseline'] > 1e-15) & \
            (cgs_to_R(drpall_geo['geo_fluxerr']) < 0.4) #& (drpall_geo['geo_sh_hght']> 500)
    sub = drpall_geo[iis0]

    
    
    x_all = sub['geo_sh_hght'].data
    y_all = cgs_to_R(sub['geo_flux'].data) 
    y_err = cgs_to_R(sub['geo_fluxerr'].data)

    # Initial guesses
    tau_initial = 3.0  # Fixed scale length
    initial_r0 = 5  # Initial guess for r0
    initial_rbgr = 2  # Initial guess for background
    initial_params = np.concatenate([[tau_initial], [initial_r0], [initial_rbgr]])

    # Fit the data
    result = least_squares(
        residuals_global_errs, 
        initial_params, 
        args=(x_all, y_all, y_err),
    )
    return result.x
    # return tau, r0, rbgr for global fit
    
def create_geocoronal_night_models(ver,save_lookup = True):
# for a given DRP version, return a lookup table of parameters for a geocoronal model, 
# taking half nights 
# returns a lookup table of the fit parameters per night

    drpall_geo = Table.read('drpall-'+ver+'_geo.fits')

    # figure out the halves of the night
    obsT = Time(drpall_geo['obstime'],format='isot', scale='utc') 
    time = obsT.jd % 1 * 24 
    time[time>0]-=24

    sunset = (time < -7.5) 
    dawn = (time >= -7.5)

    sections = [sunset, dawn]
    titles = ['start (t <-7.5h)','end (-7.5 < t) ']
    
    lookup_tables = []
    lookup_table_start = {}
    lookup_table_end = {}        

    for ii,ss in enumerate(sections):
        # no cloud issues, based on sky continuum
        #low errors on the geocoronal line fit
        iis0 =  (drpall_geo['moon_fli'] < 0.2) & (drpall_geo['geo_baseline'] < 7e-15) & (drpall_geo['geo_baseline'] > 1e-15) & \
                (cgs_to_R(drpall_geo['geo_fluxerr']) < 0.4) & \
                ss & (drpall_geo['geo_sh_hght'] < 15000) #& (t0['mjd']<60400)

        # use only data where the moon is completely down?
#        iis0 =  (drpall_geo['moon_alt'] < 0) & (drpall_geo['geo_baseline'] < 7e-15) & (drpall_geo['geo_baseline'] > 1e-15) & \
#                (cgs_to_R(drpall_geo['geo_fluxerr']) < 0.4) & \
#                ss & (drpall_geo['geo_sh_hght'] < 15000) #& (t0['mjd']<60400)

        t = drpall_geo[iis0]

        # Count the number of measurements for each unique night
        night_counts = Counter(t['mjd'])

        # Create a boolean array: True if a night has more than 3 measurements, otherwise False
        night_has_more_than_5 = np.array([night_counts[n] > 5 for n in t['mjd']])


        iis =  night_has_more_than_5 
        sub = t[iis]

        x_all = sub['geo_sh_hght'].data
        y_all = cgs_to_R(sub['geo_flux'].data) 
        y_err = cgs_to_R(sub['geo_fluxerr'].data)
        y_res = y_all * np.nan

        nights = sub['mjd'].data    

        # Step 1: Get the unique values of nights
        unique_nights = np.unique(nights)

        # Step 2: Create a mapping from unique values to integers
        night_to_index = {night: idx for idx, night in enumerate(unique_nights)}

        # Step 3: Map the original nights array to indices
        nights_indices = np.array([night_to_index[night] for night in nights])


        # Initial guesses
        tau_initial = 3.0  # Fixed scale length
        num_nights = len(np.unique(nights))
        initial_r0 = np.ones(num_nights) * 5  # Initial guess for r0
        initial_rbgr = np.zeros(num_nights)+2  # Initial guess for background
        initial_params = np.concatenate([[tau_initial], initial_r0, initial_rbgr])

        # Define bounds (all parameters must be positive)
        lower_bounds = np.concatenate([[0], initial_r0*0, initial_rbgr*0])  # all parameters >= 0
        upper_bounds = [np.inf] * len(initial_params)  # No upper limit

        # Fit the data
    #    result = least_squares(residuals, initial_params, args=(x_all, y_all, nights_indices))
        result = least_squares(
            residuals_errs, 
            initial_params, 
            args=(x_all, y_all, nights_indices,y_err),
            bounds=(lower_bounds, upper_bounds)
        )
        # Extract fitted parameters
        fitted_tau = result.x[0]
        fitted_r0 = result.x[1:num_nights+1]
        fitted_rbgr = result.x[num_nights+1:]

        if (ii == 0):
            for i,night in enumerate(unique_nights):
                if fitted_r0[i] < 50:
                    lookup_table_start[night] = {'mjd':night, 
                                                 'fitted_tau':fitted_tau, 
                                                 'fitted_r0':fitted_r0[i], 
                                                 'fitted_rbgr':fitted_rbgr[i], 
                                                 'title':titles[ii],
                                                 'time':'start'}
            lookup_tables.append(lookup_table_start)

        if (ii == 1):
            for i,night in enumerate(unique_nights):
                if fitted_r0[i] < 50:
                    lookup_table_end[night] = {'mjd':night, 
                                               'fitted_tau':fitted_tau, 
                                               'fitted_r0':fitted_r0[i], 
                                               'fitted_rbgr':fitted_rbgr[i], 
                                               'title':titles[ii],
                                               'time':'end'}
            lookup_tables.append(lookup_table_end)
        
    # some have a poor fit on rbgr, replace rbgr with the median value and refit r0
    rr = []
    for nn in lookup_tables[0]:
        if (lookup_tables[0][nn]['fitted_rbgr'] > 0.1):
            rr.append(lookup_tables[0][nn]['fitted_rbgr'])
    rbgr_start=np.median(np.array(rr))
    rr = []
    for nn in lookup_tables[1]:
        if (lookup_tables[1][nn]['fitted_rbgr'] > 0.1):
            rr.append(lookup_tables[1][nn]['fitted_rbgr'])
    rbgr_end=np.median(np.array(rr))

    for lookup_table in lookup_tables:
        for nn in lookup_table:
            fitted_tau = lookup_table[nn]['fitted_tau']
            fitted_r0 = lookup_table[nn]['fitted_r0']
            fitted_rbgr_orig = lookup_table[nn]['fitted_rbgr']

            if fitted_rbgr_orig < 0.1:
                if lookup_table[nn]['time'] == 'start':
                    lookup_table[nn]['fitted_rbgr'] = rbgr_start
                    ss = sunset
                else:
                    lookup_table[nn]['fitted_rbgr'] = rbgr_end
                    ss = dawn
                fitted_rbgr = lookup_table[nn]['fitted_rbgr'] 

                # Initial guesses
                initial_r0 =  5  # Initial guess for r0

                iis0 =  (drpall_geo['moon_alt']<0) & (drpall_geo['geo_baseline'] > 1e-15) & \
                        (cgs_to_R(drpall_geo['geo_fluxerr']) < 0.4) & \
                        ss & (drpall_geo['geo_sh_hght']< 15000) & (drpall_geo['mjd']==nn)

                t = drpall_geo[iis0]
                x_all = t['geo_sh_hght']
                y_all = cgs_to_R(t['geo_flux'].data) 
                y_err = cgs_to_R(t['geo_fluxerr'].data)


                # Fit the data
                result = least_squares(
                    residuals_errs_fixrbgr, 
                    initial_r0, 
                    args=(x_all, y_all, fitted_tau,fitted_rbgr ,y_err)
                )        
                new_fitted_r0 = result.x[0]
                lookup_table[nn]['fitted_r0'] = new_fitted_r0


    if save_lookup:
        pickle.dump( lookup_tables, open( "geocoronal_lookup_tables_"+ver+".pickle", "wb" ) )

    return lookup_tables    
    # return lookup table

def predict_geocoronal(ver, obstime, sh,force_neighbors=False, force_global=False, create_lookup_table=False):
# returns the modeled geocoronal Ha flux, given a specific DRP version, time of observation and shadow height
# ver = '1.1.0'
# obstime = '2023-11-08T02:15:23.174'  # format='isot', in UTC
# sh = 968.6                           # in km
# force_neighbors = True               # uses the lookup table fits from neighboring nights
# force_global = True                  # forces it to use only simple the global model
# create_lookup_table = True           # forces it to create the lookup table (takes a few seconds)
# 
# returns: 
# model_geocoronal            # the modeled geocoronal value in erg/s/cm^2/arcsec^2
# quality_flag                # 1 = global model, 2 = neighboring nights, 3 = model from specific night

    if force_global:
        params = create_geocoronal_global_model(ver)
        model_geocoronal = model_global(params,sh)
        quality_flag = 1
        return model_geocoronal, quality_flag
    
    if (os.path.exists('geocoronal_lookup_tables_'+ver+'.pickle')) & (not create_lookup_table):
        lookup_tables = pickle.load( open( 'geocoronal_lookup_tables_'+ver+'.pickle', "rb" ) )
    else:
        lookup_tables = create_geocoronal_night_models(ver, save_lookup = True) # saves time to save the table
    
    obsT = Time(obstime,format='isot', scale='utc') 
    tt = obsT.jd % 1 * 24 
    nn = np.floor(obsT.mjd)
    model_geocoronal = np.nan

    if sh > 15000:
#        print('shadow height too large, not constrained by the model')
        params = create_geocoronal_global_model(ver)
        model_geocoronal = model_global(params,sh)
        quality_flag = 1
        return model_geocoronal, quality_flag

    # figure out if first or 2nd half of the night
    if tt > 0:
        tt -= 24
    if (tt < -7.5): #sunset 
        lookup_table = lookup_tables[0]
    else: # dawn
        lookup_table = lookup_tables[1]

    # which night to use?    
    if (not force_neighbors) and (nn in lookup_table):    
        fitted_tau = lookup_table[nn]['fitted_tau']
        fitted_r0 = lookup_table[nn]['fitted_r0']
        fitted_rbgr = lookup_table[nn]['fitted_rbgr']
        quality_flag = 3
    else:
        quality_flag = 2
        # make this an iterative attempt to find close neighboring nights
        fitted_tau_arr=[]
        fitted_r0_arr=[]
        fitted_rbgr_arr=[]
        nn_arr = [nn+i for i in range(-2,3)]
        for nn1 in nn_arr:
            if nn1 in lookup_table:    
                fitted_tau_arr.append(lookup_table[nn1]['fitted_tau'])
                fitted_r0_arr.append(lookup_table[nn1]['fitted_r0'])
                fitted_rbgr_arr.append(lookup_table[nn1]['fitted_rbgr'])
        if (len(fitted_tau_arr)> 0):
            fitted_tau = np.median(np.array(fitted_tau_arr))
            fitted_r0 = np.median(np.array(fitted_r0_arr))
            fitted_rbgr = np.median(np.array(fitted_rbgr_arr))
        else:
            nn_arr = [nn+i for i in range(-5,6)]
            fitted_tau_arr=[]
            fitted_r0_arr=[]
            fitted_rbgr_arr=[]
            for nn1 in nn_arr:
                if nn1 in lookup_table:    
                    fitted_tau_arr.append(lookup_table[nn1]['fitted_tau'])
                    fitted_r0_arr.append(lookup_table[nn1]['fitted_r0'])
                    fitted_rbgr_arr.append(lookup_table[nn1]['fitted_rbgr'])
            if (len(fitted_tau_arr)> 0):
                fitted_tau = np.median(np.array(fitted_tau_arr))
                fitted_r0 = np.median(np.array(fitted_r0_arr))
                fitted_rbgr = np.median(np.array(fitted_rbgr_arr))
            else:
                nn_arr = [nn+i for i in range(-50,51)]
                fitted_tau_arr=[]
                fitted_r0_arr=[]
                fitted_rbgr_arr=[]
                for nn1 in nn_arr:
                    if nn1 in lookup_table:    
                        fitted_tau_arr.append(lookup_table[nn1]['fitted_tau'])
                        fitted_r0_arr.append(lookup_table[nn1]['fitted_r0'])
                        fitted_rbgr_arr.append(lookup_table[nn1]['fitted_rbgr'])
                if (len(fitted_tau_arr)> 0):
                    fitted_tau = np.median(np.array(fitted_tau_arr))
                    fitted_r0 = np.median(np.array(fitted_r0_arr))
                    fitted_rbgr = np.median(np.array(fitted_rbgr_arr))
                else:
#                    print('no data from neighboring nights')
                    params = create_geocoronal_global_model(ver)
                    model_geocoronal = model_global(params,sh)
                    quality_flag = 1
                    return model_geocoronal, quality_flag

    params = np.concatenate([[fitted_tau], [fitted_r0], [fitted_rbgr]])
    model_geocoronal = model_global(params, sh) #fitted_r0 * np.exp(-1.*(sh/1000.) / fitted_tau) + fitted_rbgr
    if model_geocoronal > 20:
#        print('unphysical modeled Ha')
        params = create_geocoronal_global_model(ver)
        model_geocoronal = model_global(params,sh)
        quality_flag = 1
        return model_geocoronal, quality_flag 
    
    return model_geocoronal, quality_flag
