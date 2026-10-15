"""Representative samples retain stress cases and add genuinely adopted fits."""
import numpy as np
from astropy.io import fits
from astropy.table import Table
from report_cases import diagnostic_cases,case_table


def test_three_additional_examples_use_adopted_z_when_r_is_survey():
    wave=np.arange(3600.,9801.,2);n=14
    lines=np.maximum(np.sin(wave/4),0)**12*1e-15
    observed=np.tile(3e-14+lines,(n,1))
    qa=Table(rows=[dict(input_row=i,channel=ch,accepted=ch=='Z',reason='insitu' if ch=='Z' else 'line_information',
                       quality_bin_counts=[3,4,5,6] if ch=='Z' else [0,0,0,1]) for i in range(n) for ch in 'RZ'])
    h=fits.HDUList([fits.PrimaryHDU(),fits.ImageHDU(wave,name='WAVE'),
        fits.BinTableHDU(Table(dict(input_row=np.arange(n))),name='LSF_FIBERS'),
        fits.BinTableHDU(Table(dict(fiberid=np.arange(n)+1)),name='SLITMAP'),
        fits.BinTableHDU(Table(dict(ratio_gaia_to_lvm_flux_plus_sky_g=np.arange(n)*.002,gaia_g_flux=np.arange(n)*1e-14)),name='GAIA'),
        fits.BinTableHDU(qa,name='LSF_QA')])
    for name,array in [('OBSFLUX',observed),('FIT_MODEL',observed+np.arange(1,n+1)[:,None]*1e-16),
        ('IVAR',np.ones_like(observed)*1e30),('FIT_MASK',np.zeros_like(observed,dtype='u1')),
        ('MASK',np.zeros_like(observed,dtype='u1')),('FIT_LINES',np.tile(lines,(n,1))),
        ('FIT_CONT',np.zeros_like(observed)),('FIT_SOLAR',np.ones_like(observed)*3e-14),
        ('MED_FIT_MASK',np.zeros(len(wave),dtype='u1'))]:h.append(fits.ImageHDU(array,name=name))
    cases=diagnostic_cases(h,16819)
    assert len(cases)==len({c['fiberid'] for c in cases})==8
    extra=[c for c in cases if c['selection']=='Additional adopted in-situ example']
    assert len(extra)==3 and all(c['qa']['Z']['accepted'] and not c['qa']['R']['accepted'] for c in extra)
    table=case_table(cases,16819,wave)
    assert 'Q1: 3, Q2: 4, Q3: 5, Q4: 6' in table and 'survey — insufficient line coverage' in table
    assert '0/4 regions have ≥2 lines → survey' in table and '4/4 regions have ≥2 lines → in situ' in table
    assert 'Q1: 5788.0–6204.5' in table
    assert 'Q1–Q4 are consecutive wavelength regions, blue → red' in table
    assert 'counts' in table.lower() and 'not acceptance thresholds' in table

    assert '<caption' in table and '[1, 0, 2, 6] passes 2/4 regions → survey LSF' in table
