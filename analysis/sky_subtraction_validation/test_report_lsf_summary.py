"""Count fibers once across channels and keep invalid channels separate."""
import numpy as np
import pytest
from astropy.io import fits
from astropy.table import Table
from report_lsf_summary import summarize_fallbacks,fallback_table


def product():
    rows=[dict(input_row=i,channel=ch,accepted=accepted,reason=reason) for i,ch,accepted,reason in
          [(0,'R',False,'line_information'),(1,'R',False,'profile_guard'),(2,'R',True,'insitu'),
           (0,'Z',False,'line_information'),(1,'Z',True,'insitu'),(2,'Z',False,'unusable_survey_fallback'),
           (3,'R',False,'line_information'),(3,'Z',False,'line_information')]]
    return fits.HDUList([fits.PrimaryHDU(),
        fits.BinTableHDU(Table(dict(telescope=['Sci','Sci','Sci','SkyE','Sci'])),name='SLITMAP'),
        fits.BinTableHDU(Table(dict(input_row=[0,1,2,3])),name='LSF_FIBERS'),
        fits.BinTableHDU(Table(rows=rows),name='LSF_QA')])


def test_summary_counts_only_processed_science_and_deduplicates_channels():
    h=product();row=summarize_fallbacks(h,123)
    assert row['processed']==3 and row['science']==4 and row['survey_either_channel']==2
    assert row['channels']['R']==dict(insitu=1,survey=2,invalid=0,reasons={'line_information':1,'profile_guard':1})
    assert row['channels']['Z']==dict(insitu=1,survey=1,invalid=1,reasons={'line_information':1})
    for data in row['channels'].values():
        assert data['survey']==sum(data['reasons'].values())
        assert data['insitu']+data['survey']+data['invalid']==row['processed']
    html=fallback_table([row])
    assert '2 (66.7%)' in html and '1 (33.3%)' in html
    assert html.count('<table>')==2 and html.index('R-channel')<html.index('Z-channel')
    assert html.count('In situ · fibers (%)')==2 and html.count('Reasons for survey · fibers')==2
    assert 'Invalid channel (not survey): 1' in html and 'primary outcomes' in html


def test_duplicate_or_missing_qa_record_is_reported():
    h=product();h['LSF_QA'].data['input_row'][0]=1
    with pytest.raises(ValueError,match='one LSF_QA record'):summarize_fallbacks(h,123)
