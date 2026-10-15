"""Fetch the three selected native CFrames using configured SAS authentication."""
from pathlib import Path
import json,os,subprocess
from astropy.io import fits
HERE=Path(__file__).resolve().parent
OUT=Path('/data/sas/sdsswork/users/u6058164/lvm_sky_subtraction_block_validation_drp132/inputs')
rows=json.loads((HERE/'corpus_selection.json').read_text())
OUT.mkdir(parents=True,exist_ok=True)
for row in rows:
    source=Path(row['input']);destination=OUT/source.name
    if not destination.exists():
        relative=source.relative_to('/data/sas').as_posix()
        url='https://data.sdss5.org/sas/'+relative
        pending=destination.with_suffix('.partial')
        try:
            subprocess.run(['wget','--netrc','--quiet','--timeout=60','--tries=3','-O',str(pending),url],check=True)
            assert fits.getheader(pending)['EXPOSURE']==row['expnum']
            os.chmod(pending,0o644);os.link(pending,destination)
        finally:
            pending.unlink(missing_ok=True)
    row['input']=str(destination)
    print('Ready',row['expnum'],flush=True)
(HERE/'corpus_selection.json').write_text(json.dumps(rows,indent=2))
