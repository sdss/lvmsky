"""Median/ML preparation (in memory or from a seed); infer a sky-only LSF for every block."""
from concurrent.futures import ProcessPoolExecutor
from copy import deepcopy
from dataclasses import replace
from pathlib import Path
import json
import logging
import multiprocessing as mp
import os
import tempfile
import time
import ctypes
import numpy as np
from astropy.io import fits
from astropy.table import Table
from skysub.sky_subtraction import decomp_adapter as adapter
from skysub.sky_subtraction import pipeline as common
from skysub.sky_subtraction.ml import protection_velocities
from skysub.sky_subtraction.settings import SETTINGS
from skysub.sky_decomp.reliability import misfit
import re
from skysub.sky_subtraction.io import read_cframe
from skysub.sky_subtraction.fiber_lsf import unfitted_qa
from skysub.sky_subtraction.product import lsf_hdus,lvmsframe_layout,publish_product,set_prediction_reliability
from skysub.sky_decomp import reliability as rel
from skysub.sky_subtraction.reconstruction import SkyMLReconstructor
from skysub.sky_subtraction.rendering import arm_sky,exposure_line_masks,fwhm_grid
from .lsf import FIT_CHANNELS,component_matrices,individual_coefficients,channel_tensors,evaluate_block,fit_block
from .science import evaluate_science

LOG=logging.getLogger(__name__)
_WORK=None


def _fit_sky(index):
    w=_WORK; label=str(w['slit']['telescope'][index]).strip()
    model=w['renderer'].model(label); state=w['templates'][label]
    flux=w['flux'][index]*common.FACTOR; ivar=w['ivar'][index]/common.FACTOR**2
    mask=w['mask'][index].copy()
    mask[w['renderer'].h['MED_FIT_MASK'].data.astype(bool)] |= 1
    coef,status=individual_coefficients(model,state,flux,ivar,mask,w['matrices'][label])
    return dict(input_row=int(index),fiberid=int(w['slit']['fiberid'][index]),
        spectrograph=int(w['slit']['spectrographid'][index]),telescope=label,coef=coef,status=status)


def _render_arm(index):
    w=_WORK; r=w['renderer']; label=str(w['slit']['telescope'][index]).strip()
    state=evaluate_block(w['templates'][label],w['blocks'],w['slit']['spectrographid'][index],w['slit']['fiberid'][index])
    return state,arm_sky(r.model(label),r.h['MED_COEF'].data[r.labels.index(label)],state).astype('f4')


def failed_arms(provenance,medians=None):
    """Unusable sky telescopes: the recorded flag, else (older preparations) the misfit and Moon-glare tests.

    ``medians`` is ``(wave, {label: median flux}, stack META row)`` for the glare test.
    """
    out=[]
    for label in common.SKY_ARMS:
        record=provenance['median_fits'][label]
        if 'failed' in record:
            failed=bool(record['failed'])
        else:  # preparations written before the flag
            match=re.search(r'R2=([-+0-9.eE]+|nan)',record.get('summary',''))
            failed=misfit(float(match.group(1))) if match else False
            if not failed and medians is not None:
                wave,flux,meta=medians
                other,=[l for l in common.SKY_ARMS if l!=label]
                excess=rel.arm_colour_excess(flux[label],flux[other],wave)
                failed=rel.moon_glare(excess,meta[label.lower()+'_moon_sep'],meta['moon_alt'],meta['moon_fli'])
        if failed:out.append(label)
    if len(out)==len(common.SKY_ARMS):raise ValueError('Both sky arms are unusable')
    return tuple(out)


def _keep_freed_memory():
    # glibc returns every large temporary to the kernel; faulting the pages back in cost 45% of a single-thread run.
    try:
        libc=ctypes.CDLL('libc.so.6')
        libc.mallopt(-4,0);libc.mallopt(-1,2**31-1);libc.mallopt(-2,256<<20)  # M_MMAP_MAX, M_TRIM_THRESHOLD, M_TOP_PAD
    except OSError:
        pass


def _mapped(function,indices,jobs):
    if jobs==1:
        yield from map(function,indices)
    else:
        with ProcessPoolExecutor(jobs,mp_context=mp.get_context('fork')) as pool:
            yield from pool.map(function,indices,chunksize=4)


def run(input_path, *, output=None, seed=None, ensemble=None, gaia_dir=None, jobs=8, slit_count=4, centered=False, keep_preparation=False):
    global _WORK
    _keep_freed_memory()
    started=time.perf_counter(); input_path=Path(input_path).resolve()
    output=Path(output).resolve() if output else input_path.with_name(input_path.name.replace('lvmCFrame-','lvmSFrame-skyml-block-').removesuffix('.gz'))
    if jobs<1 or slit_count<4:
        raise ValueError('jobs must be positive and slit_count must be >= 4')
    if output==input_path or not output.name.startswith('lvmSFrame-skyml-block-'):
        raise ValueError('Use lvmSFrame-skyml-block-* to protect the existing methods and DRP products')
    if output.exists():raise FileExistsError(output)
    output.parent.mkdir(parents=True,exist_ok=True)
    prepared=None; seed_path=Path(seed).resolve() if seed else None
    if seed_path==output:raise ValueError('Seed and output must differ')
    if seed_path is None:
        default=output.parent/('lvmSFrame-skyml-preparation-'+input_path.name.removeprefix('lvmCFrame-').removesuffix('.gz'))
        if default.exists():seed_path=default.resolve()
        else:
            if ensemble is None:raise ValueError('Supply an ensemble for median/ML preparation or an existing --seed product')
            # Median/ML/correction stage of the common pipeline, kept in memory; per-fiber LSF fits are off.
            prepared=common.run(input_path,ensemble,gaia_dir=gaia_dir,output=default,jobs=jobs,survey_only=True,render_sky=False,return_hdus=True)
            if keep_preparation:publish_product(prepared,default,input_path);seed_path=default.resolve()
    header,slit,wave,flux,ivar,mask,_=read_cframe(input_path)
    with SkyMLReconstructor(prepared if prepared is not None else seed_path,cframe=input_path) as r:
        if prepared is None and common.sha256(input_path)!=r.provenance['input_sha256']:
            raise ValueError('Preparation belongs to a different input CFrame')
        np.testing.assert_array_equal(slit['fiberid'],r.h['SLITMAP'].data['fiberid'])
        np.testing.assert_array_equal(wave,r.wave)
        templates={label:r.state('MED_',r.labels.index(label)) for label in r.labels}
        # The residual transfer follows the current settings, whatever the preparation recorded.
        r.settings=replace(r.settings,residual_max_angstrom=SETTINGS.residual_max_angstrom,
            residual_line_protection=SETTINGS.residual_line_protection)
        velocities=list(protection_velocities(*adapter.telescope_coords(header,'Sci')))
        r.lines=exposure_line_masks(wave,templates['Sci'],velocities,r.settings)
        for label in r.labels:
            for ch in 'BRZ':
                np.testing.assert_array_equal(templates[label].knot_vectors[ch],templates['Sci'].knot_vectors[ch])
        # A misfit sky arm describes neither the sky nor the LSF: its fibers stay out of the block fit.
        failed=failed_arms(r.provenance,(r.wave,{l:np.asarray(r.h['MED_FLUX'].data[r.labels.index(l)],float) for l in common.SKY_ARMS},r.h['MED_META'].data[0]))
        stale=bool(failed) and not any('failed' in r.provenance['median_fits'][l] for l in common.SKY_ARMS)
        if stale:
            LOG.warning('Preparation predates the failed-arm flag: its ML prediction used the unusable arm %s',', '.join(failed))
        if failed:
            LOG.warning('Sky arm %s is unusable (misfit or Moon glare); excluded from the block LSF and the residual transfer',', '.join(failed))
            r.arms=[a for a,label in zip(r.arms,r.arm_labels) if label not in failed]
            r.arm_labels=[label for label in r.arm_labels if label not in failed]
        selected=r.provenance['selected_fiberids']
        skyids=set().union(*(selected[label] for label in common.SKY_ARMS if label not in failed))
        indices=np.flatnonzero(np.isin(slit['fiberid'],list(skyids)))
        matrices={label:component_matrices(r.model(label),templates[label]) for label in {str(slit['telescope'][i]).strip() for i in indices}}
        _WORK=dict(renderer=r,slit=slit,flux=flux,ivar=ivar,mask=mask,templates=templates,matrices=matrices)
        LOG.info('Decomposing %s selected Sky fibers with fixed median LSFs',len(indices))
        block_started=time.perf_counter()
        fiber_rows=list(_mapped(_fit_sky,indices,jobs)); blocks={}
        coefficient_seconds=time.perf_counter()-block_started
        tensor_seconds=0.
        for spectrograph in sorted(set(slit['spectrographid'])):
            samples=[]
            tensor_started=time.perf_counter()
            for row in fiber_rows:
                if row['spectrograph']!=spectrograph:continue
                i=row['input_row']; model=r.model(row['telescope']); state=templates[row['telescope']]
                fm=mask[i].copy();fm[r.h['MED_FIT_MASK'].data.astype(bool)] |= 1
                samples.append(dict(position=float(slit['fiberid'][i]),wave=wave,flux=flux[i]*common.FACTOR,
                    ivar=ivar[i]/common.FACTOR**2,mask=fm,designs=channel_tensors(model,row['coef'],state,matrices[row['telescope']])))
            tensor_seconds+=time.perf_counter()-tensor_started
            if not samples:raise ValueError(f'No selected Sky fibers for spectrograph {spectrograph}')
            domain=np.asarray(slit['fiberid'])[slit['spectrographid']==spectrograph]
            for ch in FIT_CHANNELS:
                coefficient,knots,qa=fit_block(samples,templates['Sci'],ch,float(domain.min()),float(domain.max()),slit_count,centered)
                if qa['slit_rank']<slit_count:raise ValueError('Selected Sky fibers do not identify the slit basis')
                # Coefficients belong to the template wavelength basis; Science renders must use it too.
                blocks[(int(spectrograph),ch)]=dict(coef=coefficient,slit_knots=knots,qa=qa,
                    wave_knots=np.array(templates['Sci'].knot_vectors[ch],float),wave_degree=int(templates['Sci'].degrees[ch]))
                LOG.info('SP%s %s: %s fibers; chi2/pixel %.3f; constraints %.2g / %.2g',spectrograph,ch,len(samples),qa['chi2_per_pixel'],qa['normalization_error'],qa['centroid_error'])
        block_timing=dict(total_seconds=time.perf_counter()-block_started,coefficient_seconds=coefficient_seconds,
            tensor_seconds=tensor_seconds,normal_seconds=sum(b['qa']['normal_seconds'] for b in blocks.values()),
            solve_seconds=sum(b['qa']['solve_seconds'] for b in blocks.values()))
        LOG.info('Block LSF timing (seconds): %s',block_timing)
        science=np.flatnonzero((np.char.strip(np.asarray(slit['telescope']).astype(str))=='Sci')&(slit['fibstatus']==0))
        lsf_row={int(v):j for j,v in enumerate(r.h['LSF_FIBERS'].data['input_row'])}
        if not set(map(int,science))<=set(lsf_row):raise ValueError('Preparation lacks adopted B survey LSFs for some Science fibers')
        LOG.info('Evaluating %s Science fibers from the global blocks',len(science))
        method='sky_block_centered' if centered else 'sky_block'
        sky,states,widths=evaluate_science(r,blocks,slit,science,lsf_row)
        # Sky fibers: telescope median model with the block R/Z surface at the fiber; B is the telescope median LSF.
        arm_rows=np.flatnonzero(np.isin(np.char.strip(np.asarray(slit['telescope']).astype(str)),[l for l in common.SKY_ARMS if l not in failed])&(slit['fibstatus']==0))
        _WORK['blocks']=blocks
        arm=list(_mapped(_render_arm,arm_rows,jobs))
        for row,(state,spectrum) in zip(arm_rows,arm,strict=True):sky[row]=spectrum
        applied=np.concatenate([science,arm_rows])
        # Preserve the common compact schema and its coefficient-based diagnostic renderer.
        # A block product used as the seed brings its own block extensions; they are replaced too.
        remove={'FLUX','IVAR','MASK','LSF','SKY','LSF_FIBERS','LSF_COEF','LSF_KNOTS','LSF_META',
            'CAND_LSF_COEF','CAND_LSF_KNOTS','CAND_LSF_META','LSF_QA','FIT_LSF_HISTORY','SKYML_META',
            'BLOCK_LSF_META','BLOCK_LSF_COEF','BLOCK_SLIT_KNOTS','BLOCK_WAVE_KNOTS','SKY_FIBER_COEF'}
        product=fits.HDUList([x.copy() for x in r.h if x.name not in remove])
        for name,values in [('FLUX',flux.astype('f4')),('IVAR',ivar.astype('f4')),('MASK',mask.copy())]:
            image=fits.ImageHDU(values,name=name)
            image.header.update(r.h[name].header)
            product.append(image)
        product['FLUX'].data[applied]-=sky[applied]
        product['SLITMAP'].data['skyml_applied'][applied]=True
        lsf=np.array(r.h['LSF'].data,copy=True);lsf[science]=np.asarray(widths)
        for row,(state,_) in zip(arm_rows,arm,strict=True):lsf[row]=fwhm_grid(state,wave)
        product.append(fits.ImageHDU(lsf,header=r.h['LSF'].header,name='LSF'))
        product.append(fits.ImageHDU(sky,header=r.h['SKY'].header,name='SKY'))
        product.append(fits.BinTableHDU(Table(dict(input_row=science,fiberid=np.asarray(slit['fiberid'])[science])),name='LSF_FIBERS'))
        # Every block profile is the adopted one, so the candidate extensions repeat it.
        for prefix in ('','CAND_'):product.extend(lsf_hdus(prefix,states))
        qa=[dict(unfitted_qa(ch,method),fiberid=int(slit['fiberid'][i]),input_row=int(i)) for i in science for ch in FIT_CHANNELS]
        product.append(fits.BinTableHDU(Table(rows=qa),name='LSF_QA'))
        block_rows=[]
        for (sp,ch),block in blocks.items():
            block_rows.append(dict(spectrograph=sp,channel=ch,wave_degree=block['wave_degree'],
                slit_degree=3,offset_half_width=templates['Sci'].config['offset_half_width_angstrom'],**block['qa']))
        product.append(fits.BinTableHDU(Table(rows=block_rows),name='BLOCK_LSF_META'))
        product.append(fits.ImageHDU(np.stack([v['coef'] for v in blocks.values()]),name='BLOCK_LSF_COEF'))
        product.append(fits.ImageHDU(np.stack([v['slit_knots'] for v in blocks.values()]),name='BLOCK_SLIT_KNOTS'))
        product.append(fits.ImageHDU(np.stack([v['wave_knots'] for v in blocks.values()]),name='BLOCK_WAVE_KNOTS'))
        product.append(fits.BinTableHDU(Table(rows=fiber_rows),name='SKY_FIBER_COEF'))
        provenance=deepcopy(r.provenance)
        for label in common.SKY_ARMS:provenance['median_fits'][label]['failed']=label in failed
        if stale:  # the preparation's prediction used the misfit arm: no fallback took place
            diag=provenance.setdefault('prediction_diagnostics',{})
            bits=np.asarray(diag.get('reliability',[0]),np.int32)
            diag['reliability']=np.where(bits<0,bits,bits|rel.RELIABILITY_ARM_MISFIT).astype(np.int32)
        provenance['settings'].update(residual_max_angstrom=r.settings.residual_max_angstrom,
            residual_line_protection=r.settings.residual_line_protection)
        provenance.setdefault('residual_correction',{})['line_protection_velocities_km_s']=velocities
        provenance.update(block_timing=block_timing,input=str(input_path),method=method,preparation=str(seed_path) if seed_path else 'in-memory',preparation_sha256=common.sha256(seed_path) if seed_path else None,
            block_settings=dict(fitted_channels='RZ',blue_lsf='each Science fiber keeps its survey LSF from the preparation; B is not fitted',slit_basis=slit_count,slit_degree=3,coordinate='SLITMAP.fiberid within spectrograph',
                wavelength_solution='input, unchanged',constraints='nonnegative, unit integral'+(', zero first moment' if centered else '; no centroid constraint')+' at all wavelength/slit control pairs',
                regularization='none',fiber_coefficients='individual full-component QP, fixed telescope-median LSF and O2 temperature',
                fit_weights='native observed IVAR and input/protected-line masks; 5 A trimmed channel edges',
                science_lsf='R and Z evaluated exclusively from SkyE/SkyW blocks',
                sky_fiber_model='SkyE/SkyW rows: telescope median decomposition rendered with the block R/Z surface at the fiber and the median B LSF; '
                'FLUX is a residual of the fibers that build the LSF and the median, not an independent test'),
            block_code_sha256={p.name:common.sha256(p) for p in Path(__file__).parent.glob('*.py')},
            elapsed_seconds=time.perf_counter()-started)
        product.append(fits.ImageHDU(np.frombuffer(json.dumps(provenance,default=common._json).encode(),dtype='u1'),name='SKYML_META'))
        if 'reliability' in provenance.get('prediction_diagnostics',{}):
            set_prediction_reliability(product[0].header,int(np.asarray(provenance['prediction_diagnostics']['reliability']).ravel()[0]))
        product[0].header['SKYMETH']='skyml-block'; product[0].header['LSFCENT']=centered; product[0].header['BLOCKCH']='RZ'
        fd,temporary=tempfile.mkstemp(prefix='.skyblock-',suffix='.fits',dir=output.parent);os.close(fd)
        try:
            lvmsframe_layout(product,output,input_path).writeto(temporary,overwrite=True,checksum=True);os.chmod(temporary,0o644);os.link(temporary,output)
        finally:
            os.unlink(temporary);product.close();_WORK=None
    return output
