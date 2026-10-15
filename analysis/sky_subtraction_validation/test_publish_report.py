"""Check organized publication, relative dependencies, and repeat publication."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parent))
from publish_report import publish,MAIN


def test_publication_has_one_entry_and_preserves_nested_links(tmp_path):
    source=tmp_path/'source';dest=tmp_path/'public';source.mkdir();dest.mkdir()
    (source/MAIN).write_text('<a href="insitu_fit_1_2.html">Fiber</a><img src="map.png"><a href="summary.json">Data</a>')
    (source/'insitu_fit_1_2.html').write_text('<a href="index.html#top">Back</a><script src="plotly.min.js"></script>')
    (source/'plotly.min.js').write_text('local runtime')
    (source/'map.png').write_bytes(b'image bytes')
    (source/'summary.json').write_text('{}')
    (dest/'index.html').write_text('old main')
    (dest/'old_report.html').write_text('<img src="old.png">')
    (dest/'old.png').write_bytes(b'old bytes')
    publish(source,dest)
    assert sorted(x.name for x in dest.glob('*.html'))==[MAIN]
    assert not (dest/'index.html').exists()
    assert 'href="spectra/insitu_fit_1_2.html"' in (dest/MAIN).read_text()
    detail=(dest/'spectra/insitu_fit_1_2.html').read_text()
    assert f'../{MAIN}#top' in detail and '../assets/plotly.min.js' in detail
    assert '../figures/old.png' in (dest/'archive/old_report.html').read_text()
    assert (dest/'figures/old.png').read_bytes()==b'old bytes'
    publish(source,dest)
    assert sorted(x.name for x in dest.glob('*.html'))==[MAIN]


def test_dependency_error_keeps_original_files(tmp_path):
    import pytest
    source=tmp_path/'source';dest=tmp_path/'public';source.mkdir();dest.mkdir()
    (source/MAIN).write_text('<img src="unknown.png">')
    (dest/'index.html').write_text('original')
    with pytest.raises(ValueError,match='unresolved dependency'):publish(source,dest)
    assert (dest/'index.html').read_text()=='original'
