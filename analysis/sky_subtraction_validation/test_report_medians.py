"""Verify saved Science selection and per-fiber correction-before-median order."""
import numpy as np
from astropy.io import fits
from astropy.table import Table
from test_report_overview import products
from report_medians import science_median_data


def test_science_median_uses_saved_selection_and_actual_subtraction_stages():
    h,drp=products();wave=h['WAVE'].data;line=np.maximum(np.sin(wave/4),0)**12*1e-15
    observed=np.tile(3*line,(3,1));raw=np.tile(line,(3,1));oh=raw*.3;res=raw*.2
    h['FLUX'].data=observed-(raw+oh+res)
    for name,values in [('OBSFLUX',observed),('SKY_ML',raw),('SKY_OH',oh),('SKY_RES',res),('SKY',raw+oh+res)]:
        if name in h:h[name].data=values
        else:h.append(fits.ImageHDU(values,name=name))
    h.append(fits.BinTableHDU(Table(dict(telescope=['Sci','SkyE','Sci'],selected=[True,False,True])),name='GAIA'))
    h.append(fits.BinTableHDU(Table(dict(telescope=['Sci'])),name='MED_LABELS'))
    h.append(fits.BinTableHDU(Table(dict(family=['line','continuum'])),name='COMPONENT_NAMES'))
    h.append(fits.ImageHDU([3*line],name='MED_FLUX'))
    h.append(fits.ImageHDU([2.7*line],name='MED_MODEL'))
    h.append(fits.ImageHDU([[line,1.7*line]],name='MED_COMPONENTS'))
    h.append(fits.ImageHDU(np.ones((1,len(wave)))*1e30,name='MED_IVAR'))
    # Keep actual per-fiber correction stages separate.
    h['SKY_ML'].data[0]+=line*.2;h['SKY_OH'].data[0]-=line*.2
    result=science_median_data(h,drp)
    assert result['n_fibers']==2
    np.testing.assert_allclose(result['observed'],3*line*1e14)
    np.testing.assert_allclose(result['stages']['ML before corrections'],1.9*line*1e14)
    np.testing.assert_allclose(result['stages']['ML + OH'],1.7*line*1e14)
    np.testing.assert_allclose(result['stages']['Final SKYML'],1.5*line*1e14)
    np.testing.assert_allclose(result['stages']['DRP'],line*1e14)
    assert np.isfinite(result['sigma']).all()

    # A nonlinear three-fiber example rules out subtracting medians of sky models.
    for product in [h,drp]:product['SLITMAP'].data['telescope'][1]='Sci'
    h['GAIA'].data['telescope'][1]='Sci';h['GAIA'].data['selected'][1]=True
    h['OBSFLUX'].data=np.array([3*line,4*line,5*line])
    h['SKY_ML'].data=np.array([line,3*line,0*line])
    h['SKY_OH'].data=np.array([.3*line,.1*line,.8*line])
    h['SKY_RES'].data=np.tile(.2*line,(3,1))
    h['SKY'].data=h['SKY_ML'].data+h['SKY_OH'].data+h['SKY_RES'].data
    h['FLUX'].data=h['OBSFLUX'].data-h['SKY'].data
    h['MED_FLUX'].data=np.array([4*line])
    result=science_median_data(h,drp)
    np.testing.assert_allclose(result['stages']['ML before corrections'],2*line*1e14)
    np.testing.assert_allclose(result['stages']['ML + OH'],1.7*line*1e14)
    np.testing.assert_allclose(result['stages']['Final SKYML'],1.5*line*1e14)


def test_median_rms_uses_input_sigma_and_preserves_physical_misfit():
    from report_medians import science_median_metrics,median_metrics_table
    wave=np.arange(3600.,9800.,.5)
    data=dict(wave=wave,observed=np.full(len(wave),3.),model=np.ones(len(wave)),
              fit_mask=np.zeros(len(wave),bool),fit_ivar=np.full(len(wave),4.),
              stages={'DRP':np.zeros(len(wave)),'ML before corrections':np.zeros(len(wave)),
                      'ML + OH':np.zeros(len(wave)),'Final SKYML':np.zeros(len(wave))})
    # Mask an extreme outlier and keep residual scatter separate from input errors.
    data['observed'][0]=1e6;data['fit_mask'][0]=True
    for row in science_median_metrics(data):
        assert row['rms']==2 and row['median_sigma']==.5 and row['normalized_rms']==4
        assert all(value==0 for value in row['stages'].values())
    assert 'MED_IVAR' in median_metrics_table(data)


def test_local_line_displacement_diagnostic_recovers_known_shift():
    from audit_median_fits import local_line_probe
    wave=np.arange(9300.,9350.,.1)
    lines=10*np.exp(-.5*((wave-9325)/.8)**2)
    observed=1.2*np.interp(wave-.12,wave,lines)+3+.01*(wave-9325)
    probe=local_line_probe(wave,observed,lines,lines+3,np.ones(len(wave))*100,9325)
    np.testing.assert_allclose(probe['shift'],.12,atol=1e-5)
    np.testing.assert_allclose(probe['amplitude'],1.2,atol=1e-5)
    assert probe['rms_local']<1e-7 and probe['rms_saved']>.1 and probe['rms_local_zero']>.1
    assert local_line_probe(wave,observed,lines,lines,np.zeros(len(wave)),9325) is None


def test_oh_status_distinguishes_applied_correction_from_retained_ml():
    from report_medians import oh_status
    assert 'applied' in oh_status({'accepted':True})
    status=oh_status({'accepted':False,'reason':'invalid populations < 0'})
    assert 'ML retained' in status and '&lt; 0' in status
