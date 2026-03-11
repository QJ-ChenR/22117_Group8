"""
map_to_structure.py – Load a PDB structure and map variants to it.

Requires Biopython (``pip install biopython``).

Usage
-----
    python scripts/map_to_structure.py \\
        --structure data/structures/protein.pdb \\
        --variants  data/variants/variants.txt \\
        --chain     A
"""

import argparse
import warnings
from pathlib import Path

from parse_variants import load_variants

try:
    from Bio import BiopythonWarning
    from Bio.PDB import PDBParser
except ImportError as exc:  # pragma: no cover
    raise ImportError(
        "Biopython is required. Install it with: pip install biopython"
    ) from exc

# Standard three-letter to one-letter amino acid mapping
_THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}


def load_structure(pdb_path: str | Path, structure_id: str = "structure"):
    """Load a PDB file and return a Biopython Structure object.

    Parameters
    ----------
    pdb_path : str or Path
    structure_id : str
        Identifier for the structure (used internally by Biopython).

    Returns
    -------
    Bio.PDB.Structure.Structure
    """
    pdb_path = Path(pdb_path)
    if not pdb_path.exists():
        raise FileNotFoundError(f"PDB file not found: {pdb_path}")

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", BiopythonWarning)
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure(structure_id, str(pdb_path))

    return structure


def get_residue_aa(structure, chain_id: str, position: int) -> str | None:
    """Return the one-letter amino-acid code for a residue in a chain.

    Parameters
    ----------
    structure : Bio.PDB.Structure.Structure
    chain_id : str
        Chain identifier (e.g. "A").
    position : int
        Residue sequence number.

    Returns
    -------
    str or None
        One-letter code, or ``None`` if the residue is absent.
    """
    for model in structure:
        if chain_id not in [c.id for c in model]:
            return None
        chain = model[chain_id]
        res_id = (" ", position, " ")
        if chain.has_id(res_id):
            resname = chain[res_id].resname.upper()
            return _THREE_TO_ONE.get(resname)
    return None


def map_variants_to_structure(
    structure,
    variants: list[dict],
    chain_id: str,
) -> list[dict]:
    """Check each variant against the structure and annotate the result.

    Parameters
    ----------
    structure : Bio.PDB.Structure.Structure
    variants : list of dict
        Output of :func:`parse_variants.load_variants`.
    chain_id : str
        Chain to map against.

    Returns
    -------
    list of dict
        Each dict is the original variant dict extended with:

        * ``chain`` – chain identifier
        * ``structure_aa`` – amino acid found in structure at that position
          (``None`` if absent)
        * ``mapped`` – ``True`` if the WT amino acid matches the structure
        * ``warning`` – human-readable warning if ``mapped`` is ``False``
    """
    mapped = []
    for variant in variants:
        pos = variant["position"]
        struct_aa = get_residue_aa(structure, chain_id, pos)
        result = dict(variant)
        result["chain"] = chain_id

        if struct_aa is None:
            result["structure_aa"] = None
            result["mapped"] = False
            result["warning"] = (
                f"Residue {pos} not found in chain {chain_id}"
            )
        elif struct_aa != variant["wt_aa"]:
            result["structure_aa"] = struct_aa
            result["mapped"] = False
            result["warning"] = (
                f"WT mismatch at position {pos}: "
                f"expected {variant['wt_aa']}, found {struct_aa} in chain {chain_id}"
            )
        else:
            result["structure_aa"] = struct_aa
            result["mapped"] = True
            result["warning"] = None

        mapped.append(result)

    return mapped


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Map missense variants to a PDB structure."
    )
    parser.add_argument("--structure", required=True, metavar="PDB",
                        help="Path to the PDB structure file.")
    parser.add_argument("--variants", required=True, metavar="FILE",
                        help="Path to the variants file.")
    parser.add_argument("--chain", default="A", metavar="CHAIN",
                        help="Chain ID to map against (default: A).")
    args = parser.parse_args(argv)

    variants = load_variants(args.variants)
    structure = load_structure(args.structure)
    mapped = map_variants_to_structure(structure, variants, args.chain)

    ok = [v for v in mapped if v["mapped"]]
    warn = [v for v in mapped if not v["mapped"]]

    print(f"Chain {args.chain} — {len(ok)}/{len(mapped)} variants mapped successfully\n")
    print(f"{'Variant':<12} {'Chain':<7} {'Struct AA':<12} {'Mapped':<8} Warning")
    print("-" * 60)
    for v in mapped:
        struct_aa = v["structure_aa"] or "—"
        warning = v["warning"] or ""
        print(
            f"{v['raw']:<12} {v['chain']:<7} {struct_aa:<12} {str(v['mapped']):<8} {warning}"
        )

    if warn:
        print(f"\nWarning: {len(warn)} variant(s) could not be mapped.")

    return mapped


if __name__ == "__main__":
    main()
