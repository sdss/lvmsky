"""Run SkyML on one lvmCFrame twice: the deployed ensemble, and the same with the O I coefficients
replaced by the dedicated networks of ``train.py``.  Run from the lvmsky repo root:

    python -m skysub.experiments_oi_dedicated_mlp_v1.exposure CFRAME --ensemble PT --models OI_MODELS_PT \
        --out-dir DIR [--jobs 6]

Writes ``DIR/deployed/lvmSFrame-skyml-<exposure>.fits`` and ``DIR/oi_mlp/lvmSFrame-skyml-<exposure>.fits``.
Everything else in the run is unchanged: only the four O I values (and their seed spread) in
``ML_COEF`` differ before the per-fiber LSF fits, rendering and subtraction.  The O I product
carries ``OIMLP`` (models file) and ``OIMAP`` (coefficient -> network) in its primary header.
"""

import os

for _name in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'RAYON_NUM_THREADS'):
    os.environ[_name] = '1'

import argparse  # noqa: E402
import logging  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402
import torch  # noqa: E402
from astropy.io import fits  # noqa: E402

from skysub.sky_subtraction import pipeline  # noqa: E402

from .train import Net  # noqa: E402

# One network per emission mechanism, as in the notebook's grouping.
NETWORK = {'ATOM_Or': 'or', 'ATOM_Orc_OI0777': 'orc', 'ATOM_Orc_OI0845': 'orc', 'ATOM_Og': 'og'}
LOG = logging.getLogger(__name__)


def oi_prediction(pack, coef_near, coef_far, triplet, names):
    """Science O I coefficients (mean, seed std) from near/far median coefficients and the ensemble context."""
    oi = [names.index(n) for n in pack['oi_names']]
    have = [str(n) for n in triplet['ctx_names']]
    ctx = [np.asarray(triplet[f'ctx_{a}'], np.float64)[:, [have.index(n) for n in pack['ctx_names']]]
           for a in ('sci', 'near', 'far')]
    log_near = np.log10(np.maximum(coef_near[oi], pack['floor']))
    log_far = np.log10(np.maximum(coef_far[oi], pack['floor']))
    x = np.nan_to_num(np.hstack([log_near[None], log_far[None], *ctx]))
    keep = pack['keep']
    x = torch.as_tensor(np.clip((x[:, keep] - pack['median'][keep]) / pack['iqr'][keep], -10, 10), dtype=torch.float32)
    mean, std = {}, {}
    for variant in sorted(set(NETWORK.values())):
        spec = pack['variants'][variant]
        runs = []
        for state in spec['states']:
            net = Net(x.shape[1], len(spec['lines']))
            net.load_state_dict(state)
            net.eval()
            with torch.no_grad():
                runs.append(net(x).numpy()[0])
        for j, line in enumerate(spec['lines']):
            if NETWORK[line] == variant:
                k = list(pack['oi_names']).index(line)
                values = 10 ** (np.array(runs)[:, j] + log_near[k])
                mean[line], std[line] = values.mean(), values.std(ddof=1)
    return mean, std


def with_oi_networks(pack):
    """A ``predict_ifu`` that runs the deployed one, then overwrites the O I coefficients."""
    deployed = pipeline.predict_ifu

    def predict_ifu(ensemble_path, frame, median_fits):
        ml, diagnostics = deployed(ensemble_path, frame, median_fits)
        slot = pipeline.arm_slots(median_fits.failed_arms)
        east, west = (median_fits.results[slot[k]].coef for k in ('SkyE', 'SkyW'))
        near_is_east = bool(np.asarray(diagnostics['near_is_east'])[0])
        near, far = (east, west) if near_is_east else (west, east)
        mean, std = oi_prediction(pack, np.asarray(near, float), np.asarray(far, float),
                                  diagnostics['triplet'], ml.names)
        for line in mean:
            k = ml.names.index(line)
            LOG.info('%s: deployed %.4g -> O I network %.4g (near %.4g)', line, ml.coef[k], mean[line], near[k])
            ml.coef[k], ml.coef_std[k], ml.oh_coef[k] = mean[line], std[line], mean[line]
        return ml, diagnostics

    return predict_ifu


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('cframe')
    p.add_argument('--ensemble', required=True)
    p.add_argument('--models', required=True)
    p.add_argument('--out-dir', required=True)
    p.add_argument('--gaia-dir')
    p.add_argument('--jobs', type=int, default=1)
    a = p.parse_args()
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(message)s')
    expnum = int(fits.getheader(a.cframe)['EXPOSURE'])
    name = f'lvmSFrame-skyml-{expnum:08d}.fits'
    out = Path(a.out_dir)
    pipeline.run(a.cframe, a.ensemble, gaia_dir=a.gaia_dir, output=out / 'deployed' / name, jobs=a.jobs)
    pack = torch.load(a.models, weights_only=False)
    pipeline.predict_ifu = with_oi_networks(pack)
    path = pipeline.run(a.cframe, a.ensemble, gaia_dir=a.gaia_dir, output=out / 'oi_mlp' / name, jobs=a.jobs)
    with fits.open(path, mode='update') as h:
        h[0].header['OIMLP'] = (Path(a.models).name, 'O I networks set ML_COEF O I')
        h[0].header['OIMAP'] = ','.join(f'{k}:{v}' for k, v in NETWORK.items())


if __name__ == '__main__':
    main()
