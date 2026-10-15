"""Run isolated block products for the previous sample and selected lunations."""
from concurrent.futures import ThreadPoolExecutor,as_completed
from pathlib import Path
import json,subprocess,sys
ROOT=Path(__file__).resolve().parents[2]
OUT=Path('/data/sas/sdsswork/users/u6058164/lvm_sky_subtraction_block_validation_drp132')
HERE=Path(__file__).resolve().parent


def execute(row):
    row=dict(row);e=row['expnum'];output=OUT/'products'/f'lvmSFrame-skyml-block-{e:08d}.fits'
    seed=row.get('seed');args=[sys.executable,'-m','skysub.sky_subtraction_block',row['input'],'--output',str(output),'--jobs','8']
    if seed:args+=['--seed',seed]
    else:args+=['--ensemble','/home/ik52/sdss5/lvmsky/mlp_ensemble_split_zodi_cont-1.3.2.pt','--gaia-dir','/data/sas/sdsswork/users/u6058164/gaia_fiber_matching_1.3.2/gaia/fibers']
    if not output.exists():
        with (HERE/f'run_{e}.log').open('w') as stream:
            subprocess.run(args,cwd=ROOT,stdout=stream,stderr=subprocess.STDOUT,check=True)
    row['output']=str(output);return row


def main():
    old=json.loads((ROOT/'analysis/sky_subtraction_validation/run_manifest.json').read_text())
    current={r['expnum']:r for r in json.loads((ROOT/'analysis/sky_subtraction_validation/run_manifest_v2.json').read_text())}
    rows=[dict(expnum=r['expnum'],input=r['input'],drp=r['drp'],seed=current[r['expnum']]['output'],purpose=current[r['expnum']]['purpose']) for r in old]
    if '--corpus' in sys.argv:rows+=json.loads((HERE/'corpus_selection.json').read_text())
    complete=[]
    with ThreadPoolExecutor(3) as pool:
        futures={pool.submit(execute,row):row for row in rows}
        for future in as_completed(futures):
            complete.append(future.result());complete.sort(key=lambda r:r['expnum'])
            (HERE/'run_manifest.json').write_text(json.dumps(complete,indent=2))
            print('Completed',futures[future]['expnum'],flush=True)


if __name__=='__main__':main()
