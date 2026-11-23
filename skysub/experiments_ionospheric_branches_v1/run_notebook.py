"""Execute the isolated-branches notebook cell by cell, saving it after every cell, with a progress log.

    python run_notebook.py [--smoke | --alpha1]

SkyML on exposure 14771 starts with the new ensemble as soon as it is saved; the notebook waits for that
product before the exposure cell.  ``--smoke`` trains 2 seeds for 2 epochs into a scratch copy and skips the
comparison and exposure cells, to check that every other cell runs.  ``--alpha1`` trains the alpha = 1
variant (FIXED_IONO_ALPHA=1) into outputs/alpha1_run.ipynb and stops after the ensemble is saved.
"""

import os
import subprocess
import sys
import time
from pathlib import Path

import nbformat
from nbclient import NotebookClient

HERE = Path(__file__).resolve().parent
LVMSKY = HERE.parents[1]
NOTEBOOK = HERE / 'notebook_ionospheric_branches.ipynb'
OUT = HERE / 'outputs'
CFRAME = Path('/Users/ik52/obs/sas/sdsswork/lvm/spectro/redux/1.3.2_fix/1044XX/1044896/60392/lvmCFrame-00014771.fits')
PERSIST_CELL = '4cae66f94995475c8221f483976bb15a'


def log(msg):
    print(f'{time.strftime("%H:%M:%S")} {msg}', flush=True)


def main():
    smoke, alpha1 = '--smoke' in sys.argv, '--alpha1' in sys.argv
    os.environ['FIXED_IONO_ALPHA'] = '1' if alpha1 else '0'
    nb = nbformat.read(NOTEBOOK, as_version=4)
    target = OUT / 'smoke_run.ipynb' if smoke else OUT / 'alpha1_run.ipynb' if alpha1 else NOTEBOOK
    if alpha1:
        nb.cells = nb.cells[:1 + next(i for i, c in enumerate(nb.cells) if c.get('id') == PERSIST_CELL)]
    if smoke:
        cfg_cell = next(i for i, c in enumerate(nb.cells) if c.get('id') == 'd6e235bc')
        nb.cells.insert(cfg_cell + 1, nbformat.v4.new_code_cell(
            "train_cfg.update(n_epochs=2, ensemble_seeds=(42, 43))\n"
            "SAVED_ENSEMBLE_PATH = str(EXPERIMENT / 'outputs' / 'smoke_ensemble.pt')"))
    product = OUT / 'exposure_14771' / 'iono_branches' / 'lvmSFrame-skyml-00014771.fits'
    exposure = None
    client = NotebookClient(nb, timeout=None, kernel_name='python3', resources={'metadata': {'path': str(HERE)}})
    with client.setup_kernel():
        for i, cell in enumerate(nb.cells):
            src = cell.source
            if smoke and ('compare.predict_lines' in src or 'exposure_paths' in src or 'bootstrap_table' in src):
                log(f'cell {i}: skipped in smoke mode')
                continue
            if 'exposure_paths' in src:
                log('waiting for SkyML on exposure 14771')
                if exposure.wait():
                    raise RuntimeError(f'SkyML on 14771 failed ({exposure.returncode}); see the exposure log')
            log(f'cell {i} ({cell.cell_type}) start')
            client.execute_cell(cell, i)
            nbformat.write(nb, target)
            log(f'cell {i} done')
            if not smoke and not alpha1 and cell.get('id') == PERSIST_CELL:
                product.parent.mkdir(parents=True, exist_ok=True)
                exposure = subprocess.Popen(
                    [str(LVMSKY / '.venv' / 'bin' / 'python'), '-m', 'skysub.sky_subtraction', str(CFRAME),
                     '--ensemble', str(OUT / 'mlp_ensemble_iono_branches_knox.pt'),
                     '--output', str(product), '--jobs', '6'],
                    cwd=LVMSKY, stdout=open(OUT / 'exposure_14771_iono_branches.log', 'w'), stderr=subprocess.STDOUT)  # noqa: SIM115
                log('SkyML on exposure 14771 started')
    log(f'finished; notebook written to {target}')


if __name__ == '__main__':
    main()
