"""
parse_variants.py – Parse a list of missense variants.

Expected variant format: <WT_AA><residue_number><MUT_AA>
Examples: R80C, A123V, G45D

Usage
-----
    python scripts/parse_variants.py --variants data/variants/variants.txt
"""

import argparse
import re
import sys
from pathlib import Path


# One-letter amino-acid codes accepted by FoldX
_AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")

# Regex for a single missense variant token (e.g. R80C or r80c)
_VARIANT_RE = re.compile(r"^([A-Za-z])(\d+)([A-Za-z])$")


def parse_variant(token: str) -> dict:
    """Parse a single variant string into a structured dict.

    Parameters
    ----------
    token : str
        Variant in the form ``<WT_AA><position><MUT_AA>`` (e.g. "R80C").

    Returns
    -------
    dict with keys: ``wt_aa``, ``position`` (int), ``mut_aa``, ``raw``

    Raises
    ------
    ValueError
        If the token is not a valid missense variant format.
    """
    token = token.strip()
    match = _VARIANT_RE.match(token)
    if not match:
        raise ValueError(
            f"Invalid variant format: '{token}'. Expected format: <WT_AA><position><MUT_AA> (e.g. R80C)"
        )
    wt_aa = match.group(1).upper()
    position = int(match.group(2))
    mut_aa = match.group(3).upper()

    if wt_aa not in _AMINO_ACIDS:
        raise ValueError(f"Unknown wild-type amino acid '{wt_aa}' in variant '{token}'")
    if mut_aa not in _AMINO_ACIDS:
        raise ValueError(f"Unknown mutant amino acid '{mut_aa}' in variant '{token}'")
    if wt_aa == mut_aa:
        raise ValueError(
            f"Wild-type and mutant amino acid are identical ('{wt_aa}') in variant '{token}'"
        )

    return {"wt_aa": wt_aa, "position": position, "mut_aa": mut_aa, "raw": token.upper()}


def load_variants(filepath: str | Path) -> list[dict]:
    """Read a variants file and return a list of parsed variant dicts.

    Lines starting with ``#`` and blank lines are ignored.

    Parameters
    ----------
    filepath : str or Path
        Path to the variants file (one variant per line).

    Returns
    -------
    list of dict
    """
    filepath = Path(filepath)
    if not filepath.exists():
        raise FileNotFoundError(f"Variants file not found: {filepath}")

    variants = []
    with filepath.open() as fh:
        for lineno, line in enumerate(fh, start=1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            try:
                variants.append(parse_variant(line))
            except ValueError as exc:
                raise ValueError(f"Line {lineno}: {exc}") from exc

    return variants


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Parse a missense-variant list and print structured information."
    )
    parser.add_argument(
        "--variants",
        required=True,
        metavar="FILE",
        help="Path to the variants file (one variant per line, e.g. R80C).",
    )
    args = parser.parse_args(argv)

    variants = load_variants(args.variants)
    print(f"{'Variant':<12} {'WT AA':<8} {'Position':<10} {'MUT AA'}")
    print("-" * 40)
    for v in variants:
        print(f"{v['raw']:<12} {v['wt_aa']:<8} {v['position']:<10} {v['mut_aa']}")
    print(f"\nTotal: {len(variants)} variant(s)")
    return variants


if __name__ == "__main__":
    main()
