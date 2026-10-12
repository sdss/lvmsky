"""Recover LSF for Gaia-clean fibers in the additional selected exposures."""
from pathlib import Path
import json
import hashlib
import os
import sys
import numpy as np
import pandas as pd
from astropy.io import fits

ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))
from analysis.lsf_epoch_coverage.build_coverage import OUT,DATA
from analysis.solarblue_lsf_serious_corpus import build_corpus
from analysis.solarblue_all_fibers import run_all_fibers as production
from skysub import decompose_parallel as dp

# This process is restricted to the exact Gaia directory required for this analysis.
build_corpus.GAIA_DIRS=(DATA/'gaia/fibers',)
PHOTOMETRY=pd.read_parquet(OUT/'photometry.parquet').set_index(['expnum','fiberid'])
ORIGINAL_BUILD=production.build_exposure_input
ORIGINAL_SHARD=production.build_shard_input
ORIGINAL_FIT_ROW=dp._fit_worker_row
SUPPLEMENTAL=OUT/'additional_clean_fits'


def fit_clean_row(kind,idx,flux_row,ivar_row,fixed_lsf_state=None):
    # A fully filtered row needs an input-status record, not 30 empty fit cycles.
    if not np.any(ivar_row>0):
        schema=dp._telluric_decomposer(kind,idx,schema_only=True)
        result=schema.failed_input_result('No active pixels after Gaia/hardware preselection')
        assert result.fit_status!='Solved'
        return result,dp._reliability_flags(result,row_index=idx)
    return ORIGINAL_FIT_ROW(kind,idx,flux_row,ivar_row,fixed_lsf_state=fixed_lsf_state)


def apply_gaia_selection(path):
    with fits.open(path,mode='update',memmap=False) as h:
        meta=h['META'].data
        keys=pd.MultiIndex.from_arrays([np.asarray(meta['expnum'],dtype=np.int64),
                                      np.asarray(meta['fiberid'],dtype=np.int64)])
        phot=PHOTOMETRY.loc[keys]
        keep=(phot.gaia_verified & phot.gaia_clean_1pct & phot.fibstatus.eq(0)).to_numpy()
        assert len(keep)==len(meta) and phot.gaia_cache.map(lambda p:Path(p).resolve().is_relative_to(DATA.resolve())).all()
        h['IVAR_SCI'].data[~keep]=0
        h['MASK_SCI'].data[~keep]=1
        meta['fit_valid_fraction'][~keep]=0
        meta['gaia_ratio']=phot.gaia_ratio.to_numpy()
        meta['gaia_clean_1pct']=phot.gaia_clean_1pct.to_numpy()
        meta['gaia_status']=np.where(phot.gaia_verified,np.where(phot.gaia_clean_1pct,'clean','contaminated'),'unavailable')
        meta['gaia_n_sources']=phot.n_gaia_sources.to_numpy()
        assert np.all(h['IVAR_SCI'].data[~keep]==0)
        assert np.all(h['MASK_SCI'].data[~keep]==1)
        assert np.all(meta['fit_valid_fraction'][~keep]==0)
        h[0].header['GAIAPRE']=True;h[0].header['GAIAPVER']=2
        h[0].header['GAIACUT']=.01;h[0].header['GAIASRC']=str(DATA/'gaia/fibers')
        for extension in h:extension.add_checksum()
    return int(keep.sum())


def prepare_clean_exposure(record,destination):
    # Publish the resumable part only after Gaia preselection succeeds.
    pending=destination.with_name(destination.stem+'.gaia-unverified.fits')
    info=ORIGINAL_BUILD(record,pending)
    info['gaia_clean_input_fibers']=apply_gaia_selection(pending)
    os.replace(pending,destination)
    return info


def prepare_clean_shard(records,destination):
    info=ORIGINAL_SHARD(records,destination)
    with fits.open(destination,mode='update',memmap=True) as h:
        h[0].header['GAIAPRE']=True;h[0].header['GAIAPVER']=2
        h[0].header['GAIACUT']=.01;h[0].header['GAIASRC']=str(DATA/'gaia/fibers')
        h[0].add_checksum()
    return info


def main():
    selected=pd.read_csv(OUT/'exposure_epoch_mapping.csv')
    extra=selected[~selected.completed_in_frozen_run].copy()
    assert len(extra)==46 and extra.source_available.all()
    assert (extra.moon_alt>60).all() and (extra.moon_fli>.95).all() and (extra.sun_alt<-18).all()
    assert set(extra.expnum)<=set(PHOTOMETRY.reset_index().expnum)
    SUPPLEMENTAL.mkdir(exist_ok=True)
    manifest=SUPPLEMENTAL/'requested_exposures.csv'
    if not manifest.exists():extra.to_csv(manifest,index=False)
    provenance=SUPPLEMENTAL/'gaia_preselection_provenance.json'
    with (OUT/'photometry.parquet').open('rb') as stream:
        photometry_sha256=hashlib.file_digest(stream,'sha256').hexdigest()
    if provenance.exists():
        previous=json.loads(provenance.read_text()).get('photometry_sha256')
        assert previous in (None,photometry_sha256),'Photometry changed: use a new supplemental output directory'
        assert previous is not None or not list((SUPPLEMENTAL/'parameters').glob('*.fits')), 'Legacy fit provenance requires independent validation'
    provenance.write_text(json.dumps(dict(
        photometry_sha256=photometry_sha256,
        gaia_directory=str(DATA/'gaia/fibers'),gaia_photometry_provenance=str(OUT/'gaia_photometry_provenance.json'),
        criterion='verified Gaia / current synthetic FLUX+SKY G <0.01 and good hardware fiber',
        exposure_geometry='user-supplied manifest: Moon altitude >60 deg, illuminated fraction >0.95, Sun altitude <-18 deg',
        excluded_spectra='zero IVAR and full mask; retain physical fiber identity',
        model=production.FIT_MODEL,additional_exposures=len(extra)),indent=2))
    for cached_input in (SUPPLEMENTAL/'work').glob('solarblue-shard-*.fits'):
        if fits.getheader(cached_input,0).get('GAIAPVER')!=2:
            count=apply_gaia_selection(cached_input)
            print(f'Validated resumed Gaia-clean input: {cached_input.name}, {count} eligible spectra',flush=True)
        with fits.open(cached_input,memmap=True) as h:
            meta=h['META'].data
            keys=pd.MultiIndex.from_arrays([np.asarray(meta['expnum'],dtype=np.int64),np.asarray(meta['fiberid'],dtype=np.int64)])
            current=PHOTOMETRY.loc[keys]
            keep=(current.gaia_verified & current.gaia_clean_1pct & current.fibstatus.eq(0)).to_numpy()
            assert np.allclose(meta['gaia_ratio'],current.gaia_ratio.to_numpy(),equal_nan=True), 'Cached Gaia ratios differ'
            assert np.array_equal(meta['gaia_clean_1pct'],current.gaia_clean_1pct.to_numpy()), 'Cached eligibility differs'
            assert np.all(h['IVAR_SCI'].data[~keep]==0) and np.all(h['MASK_SCI'].data[~keep]==1), 'Cached exclusion weights differ'
    production.build_exposure_input=prepare_clean_exposure
    production.build_shard_input=prepare_clean_shard
    dp._fit_worker_row=fit_clean_row
    sys.argv=['complete_additional_fits','--output',str(SUPPLEMENTAL),'--selected-manifest',str(manifest),
              '--workers','125','--exposures-per-shard','16','--chunk-size','1']
    production.main()


if __name__=='__main__':main()
