"""Check every stored global evaluation and all native-product subtraction arithmetic."""
from pathlib import Path
import json,sys
import numpy as np
from astropy.io import fits
from astropy import units as u
from scipy.interpolate import BSpline
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))
from skysub.sky_subtraction_block.lsf import offset_centers
HERE=Path(__file__).resolve().parent


def main():
    records=[]
    for row in json.loads((HERE/'run_manifest.json').read_text()):
        with fits.open(row['output'],memmap=True) as h,fits.open(row['input'],memmap=True) as source:
            mapping=h['LSF_FIBERS'].data;indices=np.asarray(mapping['input_row'],int);slit=h['SLITMAP'].data
            np.testing.assert_array_equal(slit['fiberid'],source['SLITMAP'].data['fiberid'])
            np.testing.assert_array_equal(h['WAVE'].data,source['WAVE'].data)
            assert len(h['BLOCK_LSF_META'].data)==6 and set(h['BLOCK_LSF_META'].data['channel'])=={'R','Z'}
            np.testing.assert_array_equal(np.flatnonzero(slit['skyml_applied']),indices)
            maximum=0.;centroid=0.;normalization=0.
            for j,meta in enumerate(h['BLOCK_LSF_META'].data):
                coefficient=np.asarray(h['BLOCK_LSF_COEF'].data[j],float)
                assert np.isfinite(coefficient).all() and coefficient.min()>=-1e-8
                assert meta['slit_rank']==4 and meta['status'] in ('Solved','AlmostSolved')
                centers=offset_centers(coefficient.shape[0],meta['offset_half_width'])
                normalization=max(normalization,float(np.max(abs(coefficient.sum(axis=0)-1))))
                centroid=max(centroid,float(np.max(abs(np.einsum('k,kqs->qs',centers,coefficient)))))
                selected=slit['spectrographid'][indices]==meta['spectrograph']
                basis=BSpline.design_matrix(mapping['fiberid'][selected],h['BLOCK_SLIT_KNOTS'].data[j],meta['slit_degree']).toarray()
                expected=np.einsum('kqs,fs->fkq',coefficient,basis)
                actual=np.asarray(h['LSF_COEF'].data[selected,'BRZ'.index(str(meta['channel']))],float)
                np.testing.assert_array_equal(h['LSF_KNOTS'].data[selected,'BRZ'.index(str(meta['channel']))],
                    np.broadcast_to(h['BLOCK_WAVE_KNOTS'].data[j],(int(selected.sum()),len(h['BLOCK_WAVE_KNOTS'].data[j]))))
                maximum=max(maximum,float(np.max(abs(expected-actual))))
                np.testing.assert_allclose(actual,expected,atol=1e-13,rtol=0)
            assert normalization<1e-7 and centroid<1e-7
            seed=row.get('seed') or str(Path(row['output']).with_name(f"lvmSFrame-skyml-preparation-{row['expnum']:08d}.fits"))
            with fits.open(seed,memmap=True) as prepared:
                np.testing.assert_array_equal(prepared['LSF_FIBERS'].data['input_row'],mapping['input_row'])
                np.testing.assert_array_equal(prepared['LSF_COEF'].data[:,0],h['LSF_COEF'].data[:,0])
                np.testing.assert_array_equal(prepared['LSF_KNOTS'].data[:,0],h['LSF_KNOTS'].data[:,0])
            scale=u.Unit(source['FLUX'].header['BUNIT']).to(u.erg/u.s/u.cm**2/u.AA)
            observed=np.asarray(source['FLUX'].data[indices],float)*scale
            delivered=np.asarray(h['FLUX'].data[indices],float)
            sky=np.asarray(h['SKY'].data[indices],float)
            assert np.isfinite(sky).all()
            np.testing.assert_allclose((delivered+sky)*1e14,observed*1e14,atol=5e-5,rtol=3e-6,equal_nan=True)
            np.testing.assert_allclose(h['IVAR'].data[indices],source['IVAR'].data[indices]/scale**2,rtol=2e-7,equal_nan=True)
            np.testing.assert_array_equal(h['MASK'].data[indices],source['MASK'].data[indices])
            assert all(x.verify_checksum()==1 for x in h)
            records.append(dict(expnum=row['expnum'],science_fibers=len(indices),blocks=6,global_evaluation_max=maximum,
                normalization_max=normalization,centroid_max_A=centroid,flux_sum_max_A_units=float(np.nanmax(abs(delivered+sky-observed))),
                fits_checksums='all pass',blue_lsf='identical to preparation survey LSF for every Science fiber',ivar_and_mask='native input',sky_fibers=len(h['SKY_FIBER_COEF'].data)))
            print('Verified',row['expnum'],flush=True)
    assert len(records)==8
    (HERE/'product_validation.json').write_text(json.dumps(records,indent=2))


if __name__=='__main__':main()
