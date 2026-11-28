"""The ionosphere corpus file: complete columns, row-aligned with the bundled ionosphere table."""

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from mlp_predictor import ionosphere, iono_corpus  # noqa: E402

CORPUS = sorted((ionosphere.CONTEXT_TABLE_DIR / 'corpus').glob('iono_corpus_*.npz'))


def test_every_quantity_exists_for_each_arm():
    for arm in iono_corpus.ARMS:
        assert f'{arm}_hmf2' in iono_corpus.COLUMNS and f'{arm}_dvtec_event' in iono_corpus.COLUMNS
    assert all(len(v) == 3 for v in iono_corpus.COLUMNS.values())


@pytest.mark.skipif(not CORPUS, reason='no corpus file')
def test_corpus_file_matches_columns_and_ionosphere_rows():
    table, cols = iono_corpus.load(CORPUS[0])
    assert set(cols) == set(iono_corpus.COLUMNS)
    assert np.array_equal(table['obstime_mjd'], ionosphere._bundled_tables()['obstime_mjd'])
    assert np.isfinite(table['sci_vtec']).all() and np.isfinite(table['far_hmf2']).all()
    assert np.isnan(table['sci_sky_altimeter']).all()  # Science's own lines are never a feature
    assert table['in_training_set'].sum() > 10000


def test_match_takes_the_nearer_neighbour_within_the_tolerance():
    from mlp_predictor.iono_corpus import _match

    mjd = np.array([1.0, 2.0])
    assert list(_match(mjd, np.array([1 + 5e-8, 1 - 5e-8, 2 - 5e-8, 2 + 5e-8, 1.5, 1 + 2e-7]))) == [0, 0, 1, 1, -1, -1]
    assert list(_match(np.array([1.0]), np.array([1.0, 3.0]))) == [0, -1]


def test_build_context_skips_what_exists_and_rebuilds_what_follows(tmp_path, monkeypatch):
    made = []
    files = {k: tmp_path / f'{k}.npz' for k in ('iono', 'roti', 'thermo', 'corpus')}
    monkeypatch.setattr(iono_corpus, 'paths', lambda stem: files)

    def maker(key):
        def make(*args, **kwargs):
            made.append(key)
            np.savez(files[key], obstime_mjd=np.array([60000.7, 60001.6]))
        return make

    monkeypatch.setattr(iono_corpus.ionosphere, 'build_context_table', maker('iono'))
    monkeypatch.setattr(iono_corpus.roti, 'fetch', lambda nights, workers: made.append(('fetch', list(nights))))
    monkeypatch.setattr(iono_corpus.roti, 'build_table', maker('roti'))
    monkeypatch.setattr(iono_corpus.thermosphere, 'build_table', maker('thermo'))
    monkeypatch.setattr(iono_corpus, 'build', maker('corpus'))
    assert iono_corpus.build_context('root', 'stem') == files
    assert made == ['iono', ('fetch', [60000, 60001]), 'roti', 'thermo', 'corpus']
    made.clear()
    iono_corpus.build_context('root', 'stem')
    assert made == []  # everything exists
    files['thermo'].unlink()
    iono_corpus.build_context('root', 'stem')
    assert made == ['thermo', 'corpus']  # a rebuilt table forces what is built from it
    made.clear()
    iono_corpus.build_context('root', 'stem', force=True)
    assert [m if isinstance(m, str) else m[0] for m in made] == ['iono', 'fetch', 'roti', 'thermo', 'corpus']
