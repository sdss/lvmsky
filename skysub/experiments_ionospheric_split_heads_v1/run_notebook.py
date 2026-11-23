"""Execute the split-heads notebook cell by cell, saving it after every cell, with a progress log.

    python run_notebook.py split   [--smoke]   # canonical run: trains the split-head ensemble, all cells
    python run_notebook.py control             # parent model retrained here; cells up to the saved ensemble

``split`` starts SkyML on exposure 14771 with the new ensemble as soon as it is saved, and waits for the
control ensemble and that product before the cells that read them.  ``--smoke`` trains 2 seeds for 2 epochs
into a scratch copy, to check that every cell runs.
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
NOTEBOOK = HERE / 'notebook_ionospheric_split_heads.ipynb'
OUT = HERE / 'outputs'
CFRAME = Path('/Users/ik52/obs/sas/sdsswork/lvm/spectro/redux/1.3.2_fix/1044XX/1044896/60392/lvmCFrame-00014771.fits')
PERSIST_CELL = '4cae66f94995475c8221f483976bb15a'


def log(msg):
    print(f'{time.strftime("%H:%M:%S")} {msg}', flush=True)


def wait_for(what, ready):
    while not ready():
        time.sleep(60)
    log(f'{what} ready')


def main():
    mode, smoke = sys.argv[1], '--smoke' in sys.argv
    os.environ['SPLIT_IONOSPHERIC'] = '1' if mode == 'split' else '0'
    nb = nbformat.read(NOTEBOOK, as_version=4)
    target = NOTEBOOK if mode == 'split' and not smoke else OUT / f'{mode}{"_smoke" if smoke else ""}_run.ipynb'
    if smoke:
        cfg_cell = next(i for i, c in enumerate(nb.cells) if c.get('id') == 'd6e235bc')
        nb.cells.insert(cfg_cell + 1, nbformat.v4.new_code_cell(
            "train_cfg.update(n_epochs=2, ensemble_seeds=(42, 43))\n"
            "SAVED_ENSEMBLE_PATH = str(EXPERIMENT / 'outputs' / 'smoke_ensemble.pt')"))
    stop = len(nb.cells) if mode == 'split' else 1 + next(i for i, c in enumerate(nb.cells)
                                                          if c.get('id') == PERSIST_CELL)
    product = OUT / 'exposure_14771' / 'split_heads' / 'lvmSFrame-skyml-00014771.fits'
    exposure = None
    client = NotebookClient(nb, timeout=None, kernel_name='python3', resources={'metadata': {'path': str(HERE)}})
    with client.setup_kernel():
        for i, cell in enumerate(nb.cells[:stop]):
            src = cell.source
            if mode == 'split' and not smoke and 'compare.predict_lines' in src:
                wait_for('control ensemble', (OUT / 'mlp_ensemble_control_knox.pt').exists)
            if mode == 'split' and not smoke and 'exposure_paths' in src:
                log('waiting for SkyML on exposure 14771')
                exposure.wait()
                if exposure.returncode:
                    raise RuntimeError(f'SkyML on 14771 failed ({exposure.returncode}); see the exposure log')
            if smoke and ('compare.predict_lines' in src or 'exposure_paths' in src
                          or 'bootstrap_table' in src):
                log(f'cell {i}: skipped in smoke mode')
                continue
            log(f'cell {i} ({cell.cell_type}) start')
            client.execute_cell(cell, i)
            nbformat.write(nb, target)
            log(f'cell {i} done')
            if mode == 'split' and not smoke and cell.get('id') == PERSIST_CELL:
                product.parent.mkdir(parents=True, exist_ok=True)
                exposure = subprocess.Popen(
                    [str(LVMSKY / '.venv' / 'bin' / 'python'), '-m', 'skysub.sky_subtraction', str(CFRAME),
                     '--ensemble', str(OUT / 'mlp_ensemble_split_iono_heads_knox.pt'),
                     '--output', str(product), '--jobs', '6'],
                    cwd=LVMSKY, stdout=open(OUT / 'exposure_14771_split_heads.log', 'w'), stderr=subprocess.STDOUT)
                log('SkyML on exposure 14771 started')
    log(f'finished; notebook written to {target}')


if __name__ == '__main__':
    main()
