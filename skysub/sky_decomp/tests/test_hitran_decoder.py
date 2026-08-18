"""Regression tests for the PALACE/HITRAN OH label decoder."""

from pathlib import Path

import numpy as np
import pytest
from astropy.table import Table

from skysub.sky_decomp.fit import OH_GROUP_KEYS, decode_hitran_id


SKYSUB_DIR = Path(__file__).resolve().parents[2]
CANONICAL_OH = SKYSUB_DIR / "palace" / "PMD" / "pmd_popmodel_OH.dat"


@pytest.mark.parametrize(
    ("branch", "expected_lower"),
    [("O", 6), ("P", 5), ("Q", 4), ("R", 3), ("S", 2)],
)
def test_rotational_branches_decode_catalog_value_as_upper_state(
    branch,
    expected_lower,
):
    table = Table(
        {
            "ID": [f"OHXXX0{branch}R2104e"],
            "vi": [10],
            "Fi": [2],
            "Ni": [4],
            "pi": ["e"],
        }
    )

    decoded = decode_hitran_id(table)

    assert decoded["v_upper"][0] == 10
    assert decoded["N_upper"][0] == 4
    assert decoded["N_lower"][0] == expected_lower
    assert decoded["F_upper"][0] == 2
    assert decoded["branch_N"][0] == branch


def test_documented_palace_example_decodes_to_n_upper_four_n_lower_two():
    decoded = decode_hitran_id(Table({"ID": ["OHXXX0SR2104e"]}))

    assert decoded["N_upper"][0] == 4
    assert decoded["N_lower"][0] == 2
    assert not bool(decoded["is_main"][0])


@pytest.mark.parametrize(
    ("identifier", "message"),
    [
        ("OHXXX0SR2104", "exactly 13 characters"),
        ("OHXXX0XR2104e", "Unsupported HITRAN rotational branches"),
        ("OHXXX0SR2101e", "negative rotational level"),
    ],
)
def test_invalid_labels_are_rejected(identifier, message):
    with pytest.raises(ValueError, match=message):
        decode_hitran_id(Table({"ID": [identifier]}))


@pytest.mark.parametrize(
    ("column", "value"),
    [("vi", 9), ("Fi", 1), ("Ni", 5), ("pi", "f")],
)
def test_reference_column_mismatches_are_rejected(column, value):
    values = {"ID": ["OHXXX0SR2104e"], column: [value]}

    with pytest.raises(ValueError, match=column):
        decode_hitran_id(Table(values))


def test_complete_canonical_catalog_matches_all_reference_labels_and_group_count():
    table = Table.read(
        CANONICAL_OH,
        format="ascii.basic",
        guess=False,
        comment="#",
        fast_reader=False,
    )

    decoded = decode_hitran_id(table)

    assert len(decoded) == 22_058
    assert len(decoded.group_by(OH_GROUP_KEYS).groups) == 461
    np.testing.assert_array_equal(decoded["v_upper"], decoded["vi"])
    np.testing.assert_array_equal(decoded["F_upper"], decoded["Fi"])
    np.testing.assert_array_equal(decoded["N_upper"], decoded["Ni"])
    np.testing.assert_array_equal(decoded["parity"].astype(str), decoded["pi"].astype(str))
    assert np.min(decoded["N_upper"]) >= 0
    assert np.min(decoded["N_lower"]) >= 0
