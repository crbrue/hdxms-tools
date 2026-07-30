import pandas as pd
import pytest

from hdxms.workflow import MASS_COL, _merge_comparison_stats


def _stats(rows):
    return pd.DataFrame(
        rows,
        columns=[
            "Start", "End", "Sequence", "Charge", MASS_COL,
            "avg_pctD", "sd_pctD", "n",
        ],
    )


def test_comparison_uses_only_common_peptide_identities():
    a = _stats([
        (1, 5, "AAAAA", 2, 500.1234564, 30.0, 1.0, 3),
        (6, 10, "BBBBB", 3, 600.2, 50.0, 2.0, 3),
    ])
    b = _stats([
        (1, 5, "AAAAA", 2, 500.12345649, 10.0, 1.5, 3),
        (11, 15, "CCCCC", 2, 700.3, 25.0, 1.0, 3),
    ])

    result = _merge_comparison_stats(a, b, "A", "B", "A_vs_B", "no_be")

    assert len(result) == 1
    row = result.iloc[0]
    assert (row["Start"], row["End"], row["Sequence"], row["Charge"]) == (1, 5, "AAAAA", 2)
    assert row["pctD_diff"] == pytest.approx(20.0)


def test_comparison_subtracts_values_from_the_matching_peptide_only():
    a = _stats([
        (1, 5, "AAAAA", 2, 500.1, 30.0, 1.0, 3),
        (6, 10, "BBBBB", 2, 600.2, 80.0, 1.0, 3),
    ])
    # Reverse row order to ensure subtraction is identity-based, not position-based.
    b = _stats([
        (6, 10, "BBBBB", 2, 600.2, 20.0, 1.0, 3),
        (1, 5, "AAAAA", 2, 500.1, 25.0, 1.0, 3),
    ])

    result = _merge_comparison_stats(a, b, "A", "B", "A_vs_B", "be")
    diffs = dict(zip(result["Sequence"], result["pctD_diff"]))

    assert diffs == pytest.approx({"AAAAA": 5.0, "BBBBB": 60.0})


@pytest.mark.parametrize(
    "changed_row",
    [
        (1, 5, "AAAAX", 2, 500.1, 10.0, 1.0, 3),
        (1, 5, "AAAAA", 3, 500.1, 10.0, 1.0, 3),
        (1, 5, "AAAAA", 2, 500.2, 10.0, 1.0, 3),
    ],
)
def test_sequence_charge_or_mass_mismatch_is_not_compared(changed_row):
    a = _stats([(1, 5, "AAAAA", 2, 500.1, 30.0, 1.0, 3)])
    b = _stats([changed_row])

    with pytest.raises(ValueError, match="no common peptide identities"):
        _merge_comparison_stats(a, b, "A", "B", "A_vs_B", "no_be")


def test_duplicate_identity_is_rejected_before_subtraction():
    a = _stats([
        (1, 5, "AAAAA", 2, 500.1, 30.0, 1.0, 3),
        (1, 5, "AAAAA", 2, 500.1, 35.0, 1.0, 3),
    ])
    b = _stats([(1, 5, "AAAAA", 2, 500.1, 10.0, 1.0, 3)])

    with pytest.raises(pd.errors.MergeError):
        _merge_comparison_stats(a, b, "A", "B", "A_vs_B", "be")
