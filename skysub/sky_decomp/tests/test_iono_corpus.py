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
