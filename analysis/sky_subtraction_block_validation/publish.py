"""Publish the named block report without touching other method outputs."""
from pathlib import Path
from html.parser import HTMLParser
from urllib.parse import urlsplit,unquote
import json,os,re,shutil,tempfile
HERE=Path(__file__).resolve().parent
OUT=Path('/data/sas/sdsswork/users/u6058164/lvm_sky_subtraction_block_validation_drp132')
MAIN='00_MAIN_sky_block_report.html'
ATTR=re.compile(r'''(\b(?:href|src)\s*=\s*)(["'])(.*?)(\2)''',re.I)


class Links(HTMLParser):
    def __init__(self,s):super().__init__();self.links=[];self.feed(s)
    def handle_starttag(self,tag,attrs):self.links.extend(value for key,value in attrs if key in ('href','src') and value)


def destination(name):
    if name==MAIN:return Path(name)
    suffix=Path(name).suffix
    if name.startswith('summary_') and suffix=='.html':return Path('docs')/name
    if suffix=='.html':return Path('spectra')/name
    if suffix in ('.js','.css'):return Path('assets')/name
    if suffix in ('.png','.pdf','.svg'):return Path('figures')/name
    if suffix=='.md':return Path('docs')/name
    return Path('data')/name


def publish():
    OUT.mkdir(parents=True,exist_ok=True)
    files=[HERE/MAIN,*HERE.glob('sky_fit_*.html'),*[x for x in HERE.glob('science_*.html')],*HERE.glob('median_*.html'),*HERE.glob('sky_fibers_*.html'),HERE/'plotly.min.js',
        *HERE.glob('flux_maps_*.png'),*HERE.glob('global_lsf_maps_*.png'),*HERE.glob('lsf_comparison_*.png'),*HERE.glob('lsf_profiles_*.png'),*HERE.glob('residual_budget_*.png'),*HERE.glob('individual_flux_maps_*.png'),*HERE.glob('fiber_lsf_*.png'),*HERE.glob('fiber_decomp_*.png'),*HERE.glob('fiber_decomp_spectrum_*.html'),
        *HERE.glob('*.csv'),*HERE.glob('*.json'),HERE/'README_method.md',HERE/'reviews.md',HERE/'summary_plain_english.html']
    contract=HERE.parents[1]/'skysub/sky_subtraction/FITS_STRUCTURE.md';files.append(contract)
    names={p.name for p in files}
    for source in files:
        target=OUT/destination(source.name);target.parent.mkdir(parents=True,exist_ok=True)
        fd,pending=tempfile.mkstemp(prefix='.publish-',dir=target.parent);os.close(fd)
        try:
            if source.suffix=='.html':
                def rewrite(match):
                    url=match.group(3);parts=urlsplit(url)
                    if parts.scheme or parts.netloc or not parts.path:return match.group(0)
                    name=Path(unquote(parts.path)).name
                    resolved=destination(name) if name in names else Path(parts.path)
                    if name not in names and destination(source.name).parent!=Path('.'):
                        resolved=Path(parts.path)
                    relative=os.path.relpath(OUT/resolved,target.parent)
                    if parts.fragment:relative+='#'+parts.fragment
                    return match.group(1)+match.group(2)+relative+match.group(4)
                Path(pending).write_text(ATTR.sub(rewrite,source.read_text()))
            elif source.name=='README_method.md':
                Path(pending).write_text(source.read_text().replace('../sky_subtraction/FITS_STRUCTURE.md','FITS_STRUCTURE.md'))
            else:shutil.copyfile(source,pending)
            os.chmod(pending,0o644);os.replace(pending,target)
        finally:
            Path(pending).unlink(missing_ok=True)
    for directory in [OUT,*[x for x in OUT.rglob('*') if x.is_dir()]]:directory.chmod(0o755)
    assert sorted(p.name for p in OUT.glob('*.html'))==[MAIN]
    checked=0
    for target in [OUT/MAIN,*list((OUT/'spectra').glob('*.html')),*list((OUT/'docs').glob('*.html'))]:
        text=target.read_text();assert '<details' not in text and '<summary' not in text
        for link in Links(text).links:
            parts=urlsplit(link)
            if parts.scheme or parts.netloc or not parts.path:continue
            assert (target.parent/unquote(parts.path)).resolve().is_file(),(target,link)
            checked+=1
    print(json.dumps(dict(main=str(OUT/MAIN),html_pages=1+len(list((OUT/'spectra').glob('*.html'))),checked_local_links=checked)))


if __name__=='__main__':publish()
