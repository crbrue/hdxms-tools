"""Peptide monoisotopic-mass calculations backed by Pyteomics.

The legacy SI-table convention is preserved exactly: the reported value is the
neutral monoisotopic peptide mass plus one proton mass for each observed charge.
Although the historical column is named ``Peptide monoisotopic mass (uncharged)``,
its values are charge-state-specific protonated masses, not neutral masses or m/z.
"""
from __future__ import annotations

import re
from typing import Any

PROFORMA_MARKERS = re.compile(r"[\[\]{}<>]")


def _pyteomics():
    try:
        from pyteomics import mass, parser, proforma
    except ImportError as exc:  # pragma: no cover - exercised only in broken installs
        raise RuntimeError(
            "Pyteomics is required for peptide mass calculation. Install the "
            "package with its declared dependencies, for example: pip install -e ."
        ) from exc
    return mass, parser, proforma


def _clean_text(value: Any) -> str:
    if value is None:
        return ""
    text = str(value).strip()
    return "" if text.lower() in {"", "nan", "none"} else text


def select_mass_sequence(sequence: Any, peptide_id: Any = None) -> tuple[str, bool]:
    """Choose the sequence representation used for mass calculation.

    A ProForma-like ``Peptide ID`` is preferred when it contains explicit
    modification delimiters and its unmodified sequence agrees with ``Sequence``.
    Otherwise, the plain ``Sequence`` value is used.

    Returns
    -------
    (mass_sequence, is_proforma)
    """
    plain = _clean_text(sequence).replace(" ", "")
    if not plain:
        raise ValueError("Cannot calculate peptide mass from an empty Sequence value.")

    identifier = _clean_text(peptide_id)
    if not identifier or not PROFORMA_MARKERS.search(identifier):
        return plain, False

    _, parser, proforma = _pyteomics()
    try:
        parsed = proforma.ProForma.parse(identifier)
        stripped = parser.strip(parsed)
    except Exception as exc:
        raise ValueError(
            f"Peptide ID appears to contain modifications but is not valid ProForma: {identifier!r}. "
            "Use ProForma notation such as M[Oxidation]PEPTIDE or provide a plain Sequence."
        ) from exc

    if stripped != plain:
        raise ValueError(
            "Modified Peptide ID does not match the unmodified Sequence: "
            f"Peptide ID {identifier!r} strips to {stripped!r}, but Sequence is {plain!r}."
        )
    return identifier, True


def peptide_monoisotopic_mass(
    sequence: Any,
    charge: Any,
    peptide_id: Any = None,
) -> float:
    """Return the legacy SI-table charge-state-specific monoisotopic mass.

    For an unmodified peptide this is::

        Pyteomics neutral monoisotopic mass + z * exact proton mass

    For a modified peptide represented in ProForma, the ProForma neutral mass is
    used before applying the same charge-state proton correction.
    """
    mass, _, proforma = _pyteomics()
    try:
        z = int(charge)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Invalid peptide charge state: {charge!r}") from exc
    if z < 1:
        raise ValueError(f"Peptide charge state must be >= 1, received {z}.")

    mass_sequence, is_proforma = select_mass_sequence(sequence, peptide_id)
    if is_proforma:
        neutral_mass = float(proforma.ProForma.parse(mass_sequence).mass)
    else:
        neutral_mass = float(mass.fast_mass(mass_sequence, ion_type="M", charge=0))

    proton_mass = float(mass.nist_mass["H+"][0][0])
    return neutral_mass + z * proton_mass
