"""Check row identity, paired masking, and the direction of the quality score."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
import numpy as np
from astropy.io import fits
from astropy.table import Table
from report_overview import flux_arrays,quality


def products():
    wave=np.linspace(3600,9800,3101)
    line=np.maximum(np.sin(wave/4),0)**12
    slit=Table(dict(fiberid=[3,2,1],telescope=['Sci','SkyE','Sci'],skyml_applied=[True,False,True]))
    def frame(scale):
        return fits.HDUList([fits.PrimaryHDU(),fits.BinTableHDU(slit,name='SLITMAP'),fits.ImageHDU(wave,name='WAVE'),
            fits.ImageHDU(np.tile(line*scale*1e-15,(3,1)),name='FLUX'),fits.ImageHDU(np.ones((3,len(wave)))*1e30,name='IVAR'),
            fits.ImageHDU(np.zeros((3,len(wave)),dtype='u1'),name='MASK'),fits.ImageHDU(np.tile(line*1e-15,(3,1)),name='SKY_ML'),
            fits.ImageHDU(np.zeros(len(wave),dtype='u1'),name='MED_FIT_MASK')])
    return frame(.5),frame(1.)


def test_paired_maps_keep_all_science_rows_sorted_and_share_masks():
    h,drp=products();h['MASK'].data[0,10]=1;drp['IVAR'].data[2,20]=0
    rows,good,arrays=flux_arrays(h,drp)
    np.testing.assert_array_equal(rows,[2,0])
    assert not good[1,10] and not good[0,20]
    np.testing.assert_array_equal(np.isnan(arrays[0]),np.isnan(arrays[1]))
    np.testing.assert_allclose(arrays[1][good],arrays[0][good]*.5)


def test_quality_recovers_known_residual_reduction():
    h,drp=products()
    for score in quality(h,drp,1):
        assert score['n_fibers']==2 and score['fraction_lower']==1
        np.testing.assert_allclose(score['skyml'],score['drp']*.5)


def test_report_content_is_visible_without_expanding_sections():
    from html.parser import HTMLParser
    class VisibleReport(HTMLParser):
        def handle_starttag(self,tag,attrs):
            attrs=dict(attrs)
            assert tag not in ('details','summary')
            assert 'hidden' not in attrs
            assert attrs.get('aria-expanded') != 'false'
            assert attrs.get('role') not in ('tab','tabpanel')
    root=Path(__file__).resolve().parent
    for page in [root/'00_MAIN_sky_subtraction_report.html',*root.glob('fiber_*.html'),*root.glob('insitu_fit_*.html')]:
        VisibleReport().feed(page.read_text())
    source=(root/'report_overview.py').read_text()
    assert '<details' not in source and '<summary' not in source


def test_independent_figures_sigma_guides_and_reserved_legend_space():
    import json,re
    root=Path(__file__).resolve().parent
    decoder=json.JSONDecoder()
    pages=list(root.glob('insitu_fit_*.html'))
    assert len(pages)>=25
    for page in pages:
        source=page.read_text();figures=[]
        for match in re.finditer(r'Plotly\.newPlot\(',source):
            rest=source[match.end():].lstrip();args=[]
            for _ in range(4):
                value,end=decoder.raw_decode(rest);args.append(value)
                rest=rest[end:].lstrip().removeprefix(',').lstrip()
            figures.append(args)
        assert len(figures)==2
        for _,data,layout,config in figures:
            assert config["scrollZoom"] is False
            assert layout['legend']['yref']=='container'
            assert layout['legend']['yanchor']=='bottom'
            assert layout['margin']['b']>=180 and layout['margin']['t']>=100
            assert layout['hovermode']=='closest'
            assert layout['yaxis']['title']['text'].startswith('Flux (10⁻¹⁴')
            assert not any('axis3' in key for key in layout)
        names=[t['name'] for t in figures[0][1]]
        assert names[:5]==['Fitted sky emission lines','Fitted solar/diffuse spectrum','Nuisance continuum spline','Observed spectrum','Sky + continuum fit']
        assert 'Data − model (−5 RMS)' in names
        assert figures[0][1][0]['line']['width']<figures[0][1][3]['line']['width']
        fit=figures[0][2]
        assert any(s.get('type')=='rect' and s['y0']==-1 and s['y1']==1 for s in fit['shapes'])
        assert all(any(s.get('type')=='line' and s['y0']==level and s['y1']==level for s in fit['shapes']) for level in [-1,1])
        labels={a['text']:a for a in fit['annotations']}
        assert labels['+1σ']['yshift']>0 and labels['−1σ']['yshift']<0
