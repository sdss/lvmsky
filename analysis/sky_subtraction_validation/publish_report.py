"""Publish one named report with supporting files in logical subdirectories."""
from pathlib import Path
from html.parser import HTMLParser
from urllib.parse import urlsplit,urlunsplit,unquote
import argparse
import os
import posixpath
import re
import shutil
import tempfile

MAIN='00_MAIN_sky_subtraction_report.html'
ATTR=re.compile(r'''(\b(?:href|src)\s*=\s*)(["'])(.*?)(\2)''',re.I)


class Links(HTMLParser):
    def __init__(self,text):
        super().__init__();self.links=[];self.feed(text)
    def handle_starttag(self,tag,attrs):
        self.links.extend(value for key,value in attrs if key in ('href','src') and value)


def local_name(url):
    parsed=urlsplit(url)
    if parsed.scheme or parsed.netloc or not parsed.path or parsed.path.startswith('/'):
        return None
    name=Path(unquote(parsed.path)).name
    return MAIN if name in ('index.html','index.htm') else name


def location(name):
    if name==MAIN:return Path(name)
    if name.startswith('insitu_fit_') and name.endswith('.html'):folder='spectra'
    elif name.endswith('.html'):folder='archive'
    elif Path(name).suffix in ('.png','.svg','.pdf','.webp'):folder='figures'
    elif Path(name).suffix in ('.js','.css'):folder='assets'
    elif Path(name).suffix in ('.json','.csv','.fits','.parquet'):folder='data'
    elif name.endswith('.md'):folder='docs'
    else:folder='archive'
    return Path(folder)/name


def publish(source,destination):
    source=Path(source).resolve();destination=Path(destination).resolve()
    destination.mkdir(parents=True,exist_ok=True)
    existing=[p for p in destination.rglob('*') if p.is_file()]
    files={}
    for p in existing:
        name=MAIN if p.name in ('index.html','index.htm') else p.name
        if name in files and name!=MAIN:
            raise ValueError(f'Ambiguous published filename: {name}')
        files[name]=p
    for name in list(files):
        candidate=source/name
        if candidate.is_file():files[name]=candidate
    for p in [source/MAIN,*source.glob('insitu_fit_*.html')]:
        if p.is_file():files[p.name]=p
    rules=source.parents[1]/'PLOTTING_AND_REPORT_RULES.md'
    if rules.is_file():files[rules.name]=rules
    if MAIN not in files:raise ValueError(f'Main report is required: {MAIN}')
    pending=[name for name in files if name.endswith('.html')];checked=set()
    while pending:
        name=pending.pop()
        if name in checked:continue
        checked.add(name);p=files[name]
        for url in Links(p.read_text()).links:
            dependency=local_name(url)
            if dependency is None or dependency in files:continue
            candidates=[source/dependency,(source/unquote(urlsplit(url).path)).resolve(),(p.parent/unquote(urlsplit(url).path)).resolve()]
            candidate=next((x for x in candidates if x.is_file()),None)
            if candidate is None:raise ValueError(f'{name}: unresolved dependency {url}')
            files[dependency]=candidate
            if dependency.endswith('.html'):pending.append(dependency)
    output_paths={destination/location(name) for name in files}
    for name,p in files.items():
        target=destination/location(name)
        target.parent.mkdir(parents=True,exist_ok=True)
        target.parent.chmod(0o755)
        if name.endswith('.html'):
            def rewrite(match):
                url=match.group(3);dependency=local_name(url)
                if dependency is None:return match.group(0)
                parts=urlsplit(url)
                relative=posixpath.relpath(location(dependency).as_posix(),location(name).parent.as_posix())
                updated=urlunsplit(('', '', relative, parts.query, parts.fragment))
                return match.group(1)+match.group(2)+updated+match.group(4)
            text=ATTR.sub(rewrite,p.read_text())
            fd,temporary=tempfile.mkstemp(prefix='.report-',dir=target.parent)
            try:
                with os.fdopen(fd,'w') as stream:stream.write(text)
                os.chmod(temporary,0o644);os.replace(temporary,target)
            finally:
                if Path(temporary).exists():Path(temporary).unlink()
        elif p.resolve()!=target.resolve():
            shutil.copy2(p,target)
        target.chmod(0o644)
    for p in existing:
        if p not in output_paths:p.unlink()
    for p in destination.rglob('*.html'):
        for url in Links(p.read_text()).links:
            if local_name(url) is not None:
                assert (p.parent/unquote(urlsplit(url).path)).is_file(),(p,url)
    assert sorted(p.name for p in destination.glob('*.html'))==[MAIN]
    assert not any((destination/name).exists() for name in ('index.html','index.htm'))
    print(f'Published {len(files)} files. Main report: {destination/MAIN}')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('destination',type=Path)
    args=parser.parse_args()
    publish(Path(__file__).resolve().parent,args.destination)
