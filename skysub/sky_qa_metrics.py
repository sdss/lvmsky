import numpy as np
from astropy.io import fits, ascii
from astropy.table import Table
import astropy.units as u

def run_sky_metrics(input_file='XSframe_DR20_final.fits', save_tab=True):

    out_dir = input_file[:input_file.find('fits')-1]

    with fits.open(input_file) as hdul:
        wave=hdul['WAVE'].data
        flux=hdul['FLUX'].data
        sky=hdul['SKY'].data

    #names of sky lines and cont regions, where -999 means no cont subtraction done
    sky_list=['[OI]5577', '[OI]6300', '[OH]6865', '[Na]5891', 'O2',
               'AirglowCont', 'Bcont', 'Rcont', 'Zcont']
    lrangelist=[(5575, 5580), (6298,6302), (6861,6866), (5887, 5898), (8610, 8710), 
                (5420, 5440), (4195, 4220), (6425, 6445), (9130, 9145)]
    crangelist=[(5571, 5583), (6289, 6305), (6853, 6868), (5880, 5901), (8600, 8724), 
                (-999, -999), (-999, -999), (-999, -999), (-999, -999)]

    tab=run_RSS(wave, flux, sky, sky_list, lrangelist, crangelist)

    if save_tab:
        ascii.write(tab, f"{out_dir}.dat", overwrite=True, format='fixed_width')

def sky_line_metrics(wave, flux, sky, sky_list, lrangelist, crangelist):
    resid_med_list = [] #median flux of the sky line residual
    sky_med_list = [] #median of sky flux (before sky subtraction)
    resid_sky_list = [] #residual flux over sky flux (like Knox's frac_error)

    for j in range(len(sky_list)):
        wvl, wvc = line_windows(lrangelist[j], crangelist[j], wave)
        resid_med, sky_med, resid_sky = line_stats(flux, sky, wvl, wvc)
        resid_med_list.append(resid_med)
        sky_med_list.append(sky_med)
        resid_sky_list.append(resid_sky)
    
    return resid_med_list, sky_med_list, resid_sky_list

#converts lrange and crange into wavelengths for each
def line_windows(lrange, crange, wave):
    wave_line = (wave>=lrange[0]) & (wave<=lrange[1]) # line range selection
    #ignores cont if crange is -999, and sets to -999
    if crange[0] == -999:
        wave_cont = [-999]
    else:
        wave_cont = ((wave>=crange[0]) & (wave<lrange[0])) | ((wave>lrange[1]) & (wave<=crange[1])) # cont range selection
    return wave_line, wave_cont

#calculates the sky line/region metrics
def line_stats(flux, sky, wave_line, wave_cont):
    #subtracting the average background cont if not -999
    if wave_cont[0] == -999:
        norm_line_flux = flux[wave_line]
        norm_line_sky = sky[wave_line]
    else:
        norm_line_flux = flux[wave_line]-np.nanmedian(flux[wave_cont])
        norm_line_sky = sky[wave_line]-np.nanmedian(sky[wave_cont])
        
    resid_med = np.nanmedian(norm_line_flux)
    sky_med = np.nanmedian(norm_line_sky)
    resid_sky = np.nanmedian(norm_line_flux / norm_line_sky)
    return resid_med, sky_med, resid_sky

def run_one_expnum(wave, flux, sky, sky_list, lrangelist, crangelist):
    
    resid_med_list, sky_med_list, resid_sky_list = sky_line_metrics(wave, flux, sky, sky_list, lrangelist, crangelist)

    row = []
    for j in range(len(sky_list)):
        row.append(resid_med_list[j])
        row.append(sky_med_list[j])
        row.append(resid_sky_list[j])
    return row

def run_RSS(wave, flux, sky, sky_list, lrangelist, crangelist):

    column_names = []
    for j in range(len(sky_list)):
        column_names.append(f"{sky_list[j]}_resid_med")
        column_names.append(f"{sky_list[j]}_sky_med")
        column_names.append(f"{sky_list[j]}_resid_sky")

    metrics_tab = Table(names=column_names)

    #here loop over all expnums
    for i in range(flux.shape[0]):
        new_row = run_one_expnum(wave, flux[i], sky[i], sky_list, lrangelist, crangelist)
        metrics_tab.add_row(new_row)

    return metrics_tab


if __name__ == "__main__":
    run_sky_metrics()