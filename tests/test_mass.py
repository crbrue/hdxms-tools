import math

import pytest

pyteomics = pytest.importorskip("pyteomics")
from pyteomics import mass, proforma

from hdxms.mass import peptide_monoisotopic_mass, select_mass_sequence


def test_unmodified_mass_matches_legacy_pyteomics_convention():
    sequence = "PEPTIDE"
    proton = mass.nist_mass["H+"][0][0]
    neutral = mass.fast_mass(sequence, ion_type="M", charge=0)
    for z in (1, 2, 3, 4, 5):
        observed = peptide_monoisotopic_mass(sequence, z)
        assert math.isclose(observed, neutral + z * proton, rel_tol=0, abs_tol=1e-10)


def test_higher_charge_states_are_mass_not_mz():
    sequence = "PEPTIDE"
    m2 = peptide_monoisotopic_mass(sequence, 2)
    m3 = peptide_monoisotopic_mass(sequence, 3)
    proton = mass.nist_mass["H+"][0][0]
    assert math.isclose(m3 - m2, proton, rel_tol=0, abs_tol=1e-10)


def test_proforma_modification_from_peptide_id():
    sequence = "MPEPTIDE"
    peptide_id = "M[Oxidation]PEPTIDE"
    proton = mass.nist_mass["H+"][0][0]
    expected = proforma.ProForma.parse(peptide_id).mass + 2 * proton
    observed = peptide_monoisotopic_mass(sequence, 2, peptide_id)
    assert math.isclose(observed, expected, rel_tol=0, abs_tol=1e-10)


def test_plain_nonmodified_peptide_id_does_not_replace_sequence():
    assert select_mass_sequence("PEPTIDE", "arbitrary-database-id") == ("PEPTIDE", False)


def test_mismatched_modified_peptide_id_is_rejected():
    with pytest.raises(ValueError, match="does not match"):
        peptide_monoisotopic_mass("PEPTIDE", 2, "M[Oxidation]PEPTIDE")
