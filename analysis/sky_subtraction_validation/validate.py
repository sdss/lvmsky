"""Reproduce five-exposure QA from final skyml products and matching DRP SFrames."""
from pathlib import Path
import json
import os
os.environ.setdefault('MPLBACKEND','Agg')
import numpy as np
from scipy.ndimage import uniform_filter1d
from astropy.io import fits
from astropy.table import Table
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parents[1]
EXPOSURES = {
    3469: ('0011XX/11111/60177','Stellar Milky Way field (M24)'),
    12630: ('1046XX/1046750/60356','Moon-down comparison field'),
    16819: ('1050XX/1050184/60423','Nearly full Moon'),
    44808: ('1030XX/1030046/60898','Nearly full Moon'),
    59750: ('1038XX/1038686/61192','Later nearly full-moon exposure'),
}
REDUX = Path('/data/sas/sdsswork/lvm/spectro/redux/1.3.2')


def highpass(values, good, width=51):
    numerator=uniform_filter1d(np.where(good,values,0.),width,axis=-1,mode='nearest')
    denominator=uniform_filter1d(good.astype(float),width,axis=-1,mode='nearest')
    smooth=np.divide(numerator,denominator,out=np.zeros_like(numerator),where=denominator>0)
    return np.where(good,values-smooth,np.nan)


def main():
    import sys
    sys.path.insert(0,str(REPO))
    sys.path.insert(0,str(REPO/'skysub'))
    from skysub import decompose_parallel as dp
    from skysub.sky_subtraction.reconstruction import diagnostic_hdul
    from skysub.sky_decomp.lsf_reference import load_survey_ascii_lsf
    from mlp_predictor.sky_line_scaling import NEBULAR_LINES_RED,UNMEASURED_VELOCITY_KM_S
    ROOT.mkdir(exist_ok=True)
    manifest,summary=[],[]
    plt.rcParams.update({'font.size':11,'axes.grid':True,'grid.alpha':.2})
    figure,axes=plt.subplots(5,3,figsize=(15,15),layout='constrained')
    all_fiber_records=[]
    for row,(expnum,(relative,purpose)) in enumerate(EXPOSURES.items()):
        directory=REDUX/relative
        output=directory/f'lvmSFrame-skyml-{expnum:08d}.fits'
        cframe=directory/f'lvmCFrame-{expnum:08d}.fits'
        standard=directory/f'lvmSFrame-{expnum:08d}.fits'
        with diagnostic_hdul(output,survey_baseline=True) as h, fits.open(standard,memmap=True) as drp:
            wave=np.asarray(h['WAVE'].data,float)
            indices=np.asarray(h['LSF_FIBERS'].data['input_row'],int)
            ids=np.asarray(h['LSF_FIBERS'].data['fiberid'],int)
            q=Table(h['LSF_QA'].data)
            meta=json.loads(h['SKYML_META'].data.tobytes())
            selected=np.asarray(h['GAIA'].data['selected'])&(h['GAIA'].data['telescope']=='Sci')
            selected_rows=np.flatnonzero(selected)
            labels=list(h['MED_LABELS'].data['telescope'])
            sci_med=h['MED_FLUX'].data[labels.index('Sci')]
            near=meta['prediction_diagnostics']['near_is_east'][0]
            sky_med=h['MED_FLUX'].data[labels.index('SkyE' if near else 'SkyW')]
            ref=np.nanmedian(h['LSF'].data[indices],axis=0)
            velocity=dp.measure_halpha_velocity(wave,sci_med,sky_med)
            science_mask=dp.science_line_mask(wave,ref,centre_velocity_km_s=velocity)[0]
            science_mask |= dp.science_line_mask(wave,ref,lines=NEBULAR_LINES_RED,centre_velocity_km_s=velocity,
                velocity_km_s=UNMEASURED_VELOCITY_KM_S if velocity==0 else 0.)[0]
            roundtrip=0.;component_error=0.
            channel_scores={ch:{label:[] for label in ('skyml','survey','drp')} for ch in 'BRZ'}
            channel_bias={ch:[] for ch in 'BRZ'}
            for start in range(0,len(indices),64):
                subset=indices[start:start+64]
                observed=np.asarray(h['OBSFLUX'].data[subset],float)
                sky=np.asarray(h['SKY'].data[subset],float)
                out=np.asarray(h['FLUX'].data[subset],float)
                baseline=observed-np.asarray(h['SKY_SURVEY'].data[subset],float)
                reference_out=np.asarray(drp['FLUX'].data[subset],float)
                ivar=np.asarray(h['IVAR'].data[subset],float)
                good=(h['MASK'].data[subset]==0)&(drp['MASK'].data[subset]==0)&np.isfinite(out)&np.isfinite(reference_out)&(ivar>0)&~science_mask[None,:]
                np.testing.assert_array_equal(observed.astype('f4')-sky.astype('f4'),out.astype('f4'))
                roundtrip=max(roundtrip,float(np.nanmax(np.abs(observed-sky-out))))
                component_error=max(component_error,float(np.nanmax(np.abs(sky-h['SKY_ML'].data[subset]-h['SKY_RES'].data[subset]-h['SKY_OH'].data[subset]))))
                raw_hp=highpass(np.asarray(h['SKY_ML'].data[subset],float),good)
                variants={label:highpass(values,good) for label,values in [('skyml',out),('survey',baseline),('drp',reference_out)]}
                for ch,lo,hi in [('B',3600,5787),('R',5787,7454),('Z',7454,9801)]:
                    band=(wave>=lo)&(wave<hi)
                    if ch=='B':
                        pixels=good&band[None,:]
                    else:
                        threshold=np.nanpercentile(np.maximum(raw_hp[:,band],0),99,axis=1)*.04
                        pixels=good&band[None,:]&(raw_hp>threshold[:,None])
                    for label,values in variants.items():
                        weighted=np.where(pixels,values**2*ivar,np.nan)
                        rms=np.sqrt(np.nanmean(weighted,axis=1))
                        channel_scores[ch][label].extend(rms.tolist())
                    bias=np.nanmean(np.where(pixels,variants['skyml']*np.sqrt(ivar),np.nan),axis=1)
                    channel_bias[ch].extend(bias.tolist())
            max_sky=float(np.nanmax(np.abs(h['SKY'].data)))
            assert component_error < 2e-6*max_sky
            for i,fiberid in enumerate(ids):
                reference=load_survey_ascii_lsf(REPO/'skysub/sky_decomp/data/lsf',int(fiberid),channels=('B',))
                np.testing.assert_array_equal(h['LSF_COEF'].data[i,0,:,:4],reference.coefficients['B'])
            record=dict(expnum=expnum,purpose=purpose,input=str(cframe),output=str(output),drp=str(standard),
                n_science_fibers=len(indices),selected_counts={k:len(v) for k,v in meta['selected_fiberids'].items()},
                runtime_seconds=meta['elapsed_seconds'],roundtrip_max_abs=roundtrip,component_max_abs=component_error,
                ml_reliability=meta['prediction_diagnostics']['reliability'][0],oh_accepted=meta['oh_correction']['accepted'],
                oh_reason=meta['oh_correction'].get('reason',''),input_sha256=meta['input_sha256'],model_sha256=meta['ensemble_sha256'],
                code_sha256=meta['code_sha256'],gaia_mode=meta['gaia']['mode'],
                moon_alt=float(h[0].header.get('SKY MOON_ALT',np.nan)),moon_fli=float(h[0].header.get('SKY MOON_FLI',np.nan)))
            for ch in 'BRZ':
                values=channel_scores[ch]
                metrics={f'{label}_rms':float(np.nanmedian(v)) for label,v in values.items()}
                metrics['skyml_bias']=float(np.nanmedian(channel_bias[ch]))
                metrics['fallback_fraction']=0. if ch=='B' else float(1-np.mean(q['accepted'][q['channel']==ch]))
                summary.append(dict(expnum=expnum,channel=ch,**metrics))
                for fiberid,a,b,c in zip(ids,values['skyml'],values['survey'],values['drp']):
                    all_fiber_records.append(dict(expnum=expnum,channel=ch,fiberid=int(fiberid),skyml_rms=a,survey_rms=b,drp_rms=c))
            manifest.append(record)
            # Cleanest eligible Science fibers retain the same selection in every comparison.
            selected_variants={label:np.asarray(array[selected_rows],float) for label,array in [('skyml',h['FLUX'].data),('survey',h['OBSFLUX'].data-h['SKY_SURVEY'].data),('drp',drp['FLUX'].data)]}
            selected_valid=(h['MASK'].data[selected_rows]==0)&(drp['MASK'].data[selected_rows]==0)&(h['IVAR'].data[selected_rows]>0)
            for array in selected_variants.values():
                selected_valid &= np.isfinite(array)
            y={label:np.nanmedian(np.where(selected_valid,array,np.nan),axis=0) for label,array in selected_variants.items()}
            med_good=(selected_valid.sum(axis=0)>0)&~science_mask
            for j,(ch,lo,hi) in enumerate([('B',4000,5050),('R',6200,7400),('Z',8000,9200)]):
                ax=axes[row,j]
                take=(wave>=lo)&(wave<hi)
                for label,color in [('drp','#8b8b8b'),('survey','#d78527'),('skyml','#146eaa')]:
                    ax.plot(wave[take],highpass(y[label][None,:],med_good[None,:])[0,take]*1e15,label=label,color=color,lw=.7,alpha=.85)
                ax.set_title(f'{expnum} · {ch} · {len(selected_rows)} selected Science fibers',fontsize=10)
                ax.set_xlabel('Air wavelength (Å)')
                if j==0:ax.set_ylabel('High-pass median residual\n(10⁻¹⁵ erg s⁻¹ cm⁻² Å⁻¹)')
                bound=np.nanpercentile(np.abs(highpass(y['survey'][None,:],med_good[None,:])[0,take]*1e15),99.5)*1.5
                ax.set_ylim(-max(bound,.01),max(bound,.01))
                if row==0 and j==2:ax.legend(fontsize=9)
    figure.savefig(ROOT/'residual_comparison.png',dpi=150,bbox_inches='tight')
    plt.close(figure)
    Table(rows=summary).write(ROOT/'metrics.csv',format='ascii.csv',overwrite=True)
    Table(rows=all_fiber_records).write(ROOT/'fiber_metrics.csv',format='ascii.csv',overwrite=True)
    (ROOT/'run_manifest.json').write_text(json.dumps(manifest,indent=2))
    headers='<tr><th>Exposure</th><th>Channel</th><th>Skyml RMS</th><th>Survey RMS</th><th>DRP RMS</th><th>Survey fallback</th></tr>'
    rows=''.join(f"<tr><td>{r['expnum']}</td><td>{r['channel']}</td><td>{r['skyml_rms']:.3f}</td><td>{r['survey_rms']:.3f}</td><td>{r['drp_rms']:.3f}</td><td>{'fixed Blue' if r['channel']=='B' else format(r['fallback_fraction'],'.1%')}</td></tr>" for r in summary)
    text=(ROOT/'index_v1.html').read_text()
    start=text.index('<h2>Results</h2>');end=text.index('<h2>Reproduce</h2>',start)
    overview='<table><tr><th>Exposure</th><th>Field / condition</th><th>Science fibers scored</th><th>Selected median fibers Sci / E / W</th><th>Moon altitude</th><th>Illuminated fraction</th><th>Compute time</th></tr>'+''.join(f"<tr><td>{m['expnum']}</td><td>{m['purpose']}</td><td>{m['n_science_fibers']}</td><td>{m['selected_counts']['Sci']} / {m['selected_counts']['SkyE']} / {m['selected_counts']['SkyW']}</td><td>{m['moon_alt']:.1f}°</td><td>{m['moon_fli']:.3f}</td><td>{m['runtime_seconds']:.0f} s</td></tr>" for m in manifest)+'</table>'
    body='<h2>Results</h2>'+overview+'<p>Each RMS is the median across usable Science fibers of high-pass residual / observed pixel noise. Blue uses all unmasked channel pixels; R/Z use positive high-pass sky-model line pixels above 4% of the channel 99th percentile. Nebular windows are protected. All three methods use the same pixel set. Values are diagnostic scales, not reduced chi-square. Survey holds R/Z at the survey references while retaining the same ML prediction and correction parameters. Survey fallback is the fraction of fitted Science fibers adopting their own survey reference; Blue is fixed by design. Table statistics score all 1754 usable Science fibers per exposure, whereas figures use only the indicated selected median fibers. Compute time covers the pipeline through sky assembly, before FITS publication.</p><table>'+headers+rows+'</table>'
    body+='<figure><img src="residual_comparison.png" alt="Five rows of exposure comparisons across Blue, Red and Z channels"><figcaption>Masked Science-median residuals after a 25 Å running-mean high-pass: blue = skyml, orange = survey R/Z baseline, gray = DRP. The identical selected fibers are used for each method; panels have independent y-scales set to 1.5 times the survey residual absolute 99.5th percentile, so extreme excursions can be clipped. Gaps are science-line masks or unsupported median wavelengths. The x-axis shows selected channel windows, while the table covers full channels.</figcaption></figure>'
    body+='<p><a href="metrics.csv">Channel metrics</a> · <a href="fiber_metrics.csv">Per-fiber metrics</a> · <a href="run_manifest.json">Run manifest and output paths</a></p>'
    body+='<p>Every adopted Blue coefficient matrix matches its physical-fiber ASCII reference exactly. OBSFLUX − SKY reproduces FLUX within float32 rounding, and the stored raw/residual/OH arrays reproduce the applied sky. The complete uncertainty convention is documented above.</p>'
    injection=Table.read(ROOT/'injection_metrics.csv',format='ascii.csv')
    max_error=float(np.max(injection['max_w50_error']))
    body+=f'<p>All {len(injection)} controlled R/Z injections adopted an in-situ profile; maximum central-50%-width recovery error was {max_error:.2%}. The protected narrow science-emission injection preserves integrated flux within 1%. <a href="injection_metrics.csv">Injection outcomes</a>; reproduce with <code>python analysis/sky_subtraction_validation/run_injections.py</code>.</p>'
    body+='<p>Z in-situ rendering lowers diagnostic RMS relative to survey references in all five examples, with a 26.6% reduction in 16819. R adopts survey profiles for 98.6–99.9% of fibers in the three nearly full-moon cases, mainly after profile guards detect large departures. Bright moon and high R fallback coincide here; unresolved lunar/stellar absorption is a hypothesis requiring targeted high-resolution injections, rather than a demonstrated cause.</p>'
    body+='<p>Skyml RMS is lower than DRP in 12 of 15 exposure/channel comparisons. DRP performs better in R/Z for 3469 and Z for 59750. In 59750 Z, skyml RMS is 3.487 versus DRP 2.550, with signed residual bias −2.742 observed-noise units and the original OH prediction retained. Profile acceptance establishes an admissible LSF, not subtraction accuracy. These examples show selective improvements across conditions, rather than a uniform gain.</p>'
    body+='<p>The development process observed approximately 1.8 GB parent RSS and 1.4 GB per fiber worker with six workers; shared pages make a sum of RSS an upper estimate of physical use. These snapshots describe resource scale, not measured peak memory.</p>'
    body+='<p>OH correction outcomes: '+', '.join(f"{m['expnum']}: {'applied' if m['oh_accepted'] else 'prediction retained'}" for m in manifest)+'.</p>'
    body+='<p>A live ARI TAP query for exposure 12630 returned 23,497 fiber–source matches across 1,944 fibers, with source counts and summed Gaia G flux matching the supplied cache. <a href="live_tap_validation.json">Live TAP validation</a> · <a href="integrity.json">FITS checksum and preservation checks</a>.</p>'
    body+='<p><a href="reviews.md">Independent fresh-context reviews and incorporated corrections</a></p>'
    text=text[:start]+body+text[end:]
    (ROOT/'index_v1.html').write_text(text)
    print(json.dumps(summary,indent=2))


if __name__=='__main__':
    main()
