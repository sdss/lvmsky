"""Reuse median/ML preparation; infer a sky-only LSF for every block."""
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
from skysub.sky_subtraction.rendering import exposure_line_masks
from skysub.sky_subtraction.settings import SETTINGS
from skysub.sky_subtraction.io import read_cframe
from skysub.sky_subtraction.reconstruction import SkyMLReconstructor
from skysub.sky_decomp.result_io import build_lsf_hdus
from .lsf import FIT_CHANNELS,component_matrices,individual_coefficients,channel_tensors,fit_block
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


def run(input_path, *, output=None, seed=None, ensemble=None, gaia_dir=None, jobs=8, slit_count=4, centered=False):
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
    if seed is None:
        seed=output.parent/('lvmSFrame-skyml-preparation-'+input_path.name.removeprefix('lvmCFrame-').removesuffix('.gz'))
        if not Path(seed).exists():
            if ensemble is None:raise ValueError('Supply an ensemble for median/ML preparation or an existing --seed product')
            common.run(input_path,ensemble,gaia_dir=gaia_dir,output=seed,jobs=jobs,survey_only=True,render_sky=False)
    seed=Path(seed).resolve()
    if seed==output:raise ValueError('Seed and output must differ')
    header,slit,wave,flux,ivar,mask,_=read_cframe(input_path)
    with SkyMLReconstructor(seed,cframe=input_path) as r:
        if common.sha256(input_path)!=r.provenance['input_sha256']:
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
        selected=r.provenance['selected_fiberids']
        skyids=set(selected['SkyE'])|set(selected['SkyW'])
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
        product['FLUX'].data[science]-=sky[science]
        lsf=np.array(r.h['LSF'].data,copy=True);lsf[science]=np.asarray(widths)
        product.append(fits.ImageHDU(lsf,header=r.h['LSF'].header,name='LSF'))
        product.append(fits.ImageHDU(sky,header=r.h['SKY'].header,name='SKY'))
        product.append(fits.BinTableHDU(Table(dict(input_row=science,fiberid=np.asarray(slit['fiberid'])[science])),name='LSF_FIBERS'))
        for extension in build_lsf_hdus(states):product.append(extension)
        qa=[dict(input_row=int(i),fiberid=int(slit['fiberid'][i]),channel=ch,accepted=True,
            reason='survey_fixed' if ch=='B' else method) for i in science for ch in 'BRZ']
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
        provenance['settings'].update(residual_max_angstrom=r.settings.residual_max_angstrom,
            residual_line_protection=r.settings.residual_line_protection)
        provenance.setdefault('residual_correction',{})['line_protection_velocities_km_s']=velocities
        provenance.update(block_timing=block_timing,input=str(input_path),method=method,preparation=str(seed),preparation_sha256=common.sha256(seed),
            block_settings=dict(fitted_channels='RZ',blue_lsf='each Science fiber keeps its survey LSF from the preparation; B is not fitted',slit_basis=slit_count,slit_degree=3,coordinate='SLITMAP.fiberid within spectrograph',
                wavelength_solution='input, unchanged',constraints='nonnegative, unit integral'+(', zero first moment' if centered else '; no centroid constraint')+' at all wavelength/slit control pairs',
                regularization='none',fiber_coefficients='individual full-component QP, fixed telescope-median LSF and O2 temperature',
                fit_weights='native observed IVAR and input/protected-line masks; 5 A trimmed channel edges',
                science_lsf='R and Z evaluated exclusively from SkyE/SkyW blocks'),
            block_code_sha256={p.name:common.sha256(p) for p in Path(__file__).parent.glob('*.py')},
            elapsed_seconds=time.perf_counter()-started)
        product.append(fits.ImageHDU(np.frombuffer(json.dumps(provenance,default=common._json).encode(),dtype='u1'),name='SKYML_META'))
        product[0].header['SKYMETH']='skyml-block'; product[0].header['LSFCENT']=centered; product[0].header['BLOCKCH']='RZ'
        fd,temporary=tempfile.mkstemp(prefix='.skyblock-',suffix='.fits',dir=output.parent);os.close(fd)
        try:
            product.writeto(temporary,overwrite=True,checksum=True);os.chmod(temporary,0o644);os.link(temporary,output)
        finally:
            os.unlink(temporary);product.close();_WORK=None
    return output
