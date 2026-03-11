"""
parse_foldx_output.py – Parse FoldX BuildModel and AnalyseComplex output files.

FoldX generates several output files; the ones this module reads are:

BuildModel
~~~~~~~~~~
* ``Average_<pdb>.fxout``   – averaged ΔΔG values across stochastic runs.
  Columns (tab-separated): Pdb, total energy, …
  The first data row is the WT, subsequent rows are mutants.

* ``Dif_BuildModel_<pdb>.fxout`` – ΔΔG (mutant − WT) per variant.
  Columns (tab-separated): Pdb, total energy (= ΔΔG), …

AnalyseComplex
~~~~~~~~~~~~~~
* ``Summary_AC_<pdb>.fxout`` – binding free energies.
  Contains lines with interaction energies for each chain pair.

Usage
-----
    python scripts/parse_foldx_output.py \\
        --output-dir data/results/foldx_outputs \\
        --pdb-name   protein
"""

import argparse
import re
from pathlib import Path


# ─────────────────────────────────────────────────────────────────────────────
# Internal helpers
# ─────────────────────────────────────────────────────────────────────────────

def _find_file(directory: Path, pattern: str) -> Path | None:
    """Return the first file matching *pattern* in *directory*, or ``None``."""
    matches = list(directory.glob(pattern))
    return matches[0] if matches else None


def _parse_tsv_block(lines: list[str]) -> list[dict]:
    """Parse a tab-separated block (header + data rows) into list of dicts."""
    data_lines = [l for l in lines if l.strip() and not l.startswith("//")]
    if len(data_lines) < 2:
        return []
    header = [h.strip() for h in data_lines[0].split("\t")]
    rows = []
    for line in data_lines[1:]:
        parts = [p.strip() for p in line.split("\t")]
        row = dict(zip(header, parts))
        rows.append(row)
    return rows


# ─────────────────────────────────────────────────────────────────────────────
# BuildModel output parsers
# ─────────────────────────────────────────────────────────────────────────────

def parse_dif_buildmodel(filepath: str | Path) -> list[dict]:
    """Parse a ``Dif_BuildModel_*.fxout`` file.

    Each data row corresponds to one mutation run.  The ``total energy``
    column gives ΔΔG folding (mutant energy minus WT energy), in kcal/mol.

    Parameters
    ----------
    filepath : str or Path

    Returns
    -------
    list of dict
        Keys include ``"Pdb"`` and ``"total energy"`` (ΔΔG folding).
    """
    filepath = Path(filepath)
    if not filepath.exists():
        raise FileNotFoundError(f"FoldX output file not found: {filepath}")

    with filepath.open() as fh:
        lines = fh.readlines()

    return _parse_tsv_block(lines)


def parse_average_buildmodel(filepath: str | Path) -> list[dict]:
    """Parse an ``Average_*.fxout`` file.

    Parameters
    ----------
    filepath : str or Path

    Returns
    -------
    list of dict
        All rows (WT first, then mutants).
    """
    filepath = Path(filepath)
    if not filepath.exists():
        raise FileNotFoundError(f"FoldX output file not found: {filepath}")

    with filepath.open() as fh:
        lines = fh.readlines()

    return _parse_tsv_block(lines)


def extract_ddg_folding(dif_rows: list[dict]) -> list[dict]:
    """Extract ΔΔG folding values from parsed ``Dif_BuildModel`` rows.

    Parameters
    ----------
    dif_rows : list of dict
        As returned by :func:`parse_dif_buildmodel`.

    Returns
    -------
    list of dict
        Each dict has keys ``"mutation"`` and ``"ddG_folding"`` (float, kcal/mol).
    """
    results = []
    for row in dif_rows:
        pdb_name = row.get("Pdb", "")
        # Extract mutation code from PDB name, e.g. "AR80C_1.pdb" → "AR80C"
        match = re.match(r"([A-Z][A-Z]\d+[A-Z])_?\d*\.pdb", pdb_name, re.IGNORECASE)
        mutation = match.group(1) if match else pdb_name.replace(".pdb", "")
        try:
            ddg = float(row.get("total energy", "nan"))
        except (ValueError, TypeError):
            ddg = float("nan")
        results.append({"mutation": mutation, "ddG_folding": ddg})
    return results


# ─────────────────────────────────────────────────────────────────────────────
# AnalyseComplex output parsers
# ─────────────────────────────────────────────────────────────────────────────

def parse_summary_ac(filepath: str | Path) -> list[dict]:
    """Parse a ``Summary_AC_*.fxout`` file.

    The file contains interaction energies for each PDB and chain pair.

    Parameters
    ----------
    filepath : str or Path

    Returns
    -------
    list of dict
        Keys include ``"Pdb"``, ``"Chain1"``, ``"Chain2"``, and
        ``"Interaction Energy"`` (ΔG binding, kcal/mol).
    """
    filepath = Path(filepath)
    if not filepath.exists():
        raise FileNotFoundError(f"FoldX AnalyseComplex output file not found: {filepath}")

    with filepath.open() as fh:
        lines = fh.readlines()

    return _parse_tsv_block(lines)


def compute_ddg_binding(
    wt_ac_rows: list[dict],
    mut_ac_rows: list[dict],
    interaction_col: str = "Interaction Energy",
) -> list[dict]:
    """Compute ΔΔG binding = ΔG(mutant) – ΔG(WT).

    Parameters
    ----------
    wt_ac_rows : list of dict
        AnalyseComplex rows for the wild-type structure.
    mut_ac_rows : list of dict
        AnalyseComplex rows for each mutant structure.
    interaction_col : str
        Column name for interaction energy (default: ``"Interaction Energy"``).

    Returns
    -------
    list of dict
        Each dict has keys ``"mutation"`` and ``"ddG_binding"`` (float, kcal/mol).
    """
    if not wt_ac_rows:
        return []

    try:
        wt_dg = float(wt_ac_rows[0].get(interaction_col, "nan"))
    except (ValueError, TypeError):
        wt_dg = float("nan")

    results = []
    for row in mut_ac_rows:
        pdb_name = row.get("Pdb", "")
        match = re.match(r"([A-Z][A-Z]\d+[A-Z])_?\d*\.pdb", pdb_name, re.IGNORECASE)
        mutation = match.group(1) if match else pdb_name.replace(".pdb", "")
        try:
            mut_dg = float(row.get(interaction_col, "nan"))
        except (ValueError, TypeError):
            mut_dg = float("nan")
        results.append({"mutation": mutation, "ddG_binding": mut_dg - wt_dg})

    return results


# ─────────────────────────────────────────────────────────────────────────────
# Convenience loader
# ─────────────────────────────────────────────────────────────────────────────

def load_foldx_results(
    output_dir: str | Path,
    pdb_name: str,
    is_complex: bool = False,
) -> dict:
    """Load all relevant FoldX output files from *output_dir*.

    Parameters
    ----------
    output_dir : str or Path
    pdb_name : str
        Base name of the PDB file without extension (e.g. ``"protein"``).
    is_complex : bool
        Whether to also look for AnalyseComplex results.

    Returns
    -------
    dict with keys:
        * ``"ddg_folding"`` – list of ``{mutation, ddG_folding}``
        * ``"ddg_binding"`` – list of ``{mutation, ddG_binding}`` (or ``[]``)
    """
    output_dir = Path(output_dir)

    dif_file = _find_file(output_dir, f"Dif_BuildModel_{pdb_name}.fxout")
    if dif_file is None:
        dif_file = _find_file(output_dir, f"Dif_BuildModel_{pdb_name}*.fxout")
    if dif_file is None:
        raise FileNotFoundError(
            f"No Dif_BuildModel file found for '{pdb_name}' in {output_dir}"
        )

    dif_rows = parse_dif_buildmodel(dif_file)
    ddg_folding = extract_ddg_folding(dif_rows)

    ddg_binding = []
    if is_complex:
        wt_ac_file = _find_file(output_dir, f"Summary_AC_{pdb_name}_WT.fxout")
        mut_ac_file = _find_file(output_dir, f"Summary_AC_{pdb_name}*.fxout")
        if wt_ac_file and mut_ac_file:
            wt_rows = parse_summary_ac(wt_ac_file)
            mut_rows = parse_summary_ac(mut_ac_file)
            ddg_binding = compute_ddg_binding(wt_rows, mut_rows)

    return {"ddg_folding": ddg_folding, "ddg_binding": ddg_binding}


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Parse FoldX output files and print ΔΔG values."
    )
    parser.add_argument("--output-dir", required=True, metavar="DIR",
                        help="Directory containing FoldX output files.")
    parser.add_argument("--pdb-name", required=True, metavar="NAME",
                        help="PDB base name without extension (e.g. 'protein').")
    parser.add_argument("--complex", action="store_true",
                        help="Also parse AnalyseComplex output for ΔΔG binding.")
    args = parser.parse_args(argv)

    results = load_foldx_results(args.output_dir, args.pdb_name, args.complex)

    print("ΔΔG Folding:")
    for r in results["ddg_folding"]:
        print(f"  {r['mutation']}: {r['ddG_folding']:.2f} kcal/mol")

    if results["ddg_binding"]:
        print("\nΔΔG Binding:")
        for r in results["ddg_binding"]:
            print(f"  {r['mutation']}: {r['ddG_binding']:.2f} kcal/mol")

    return results


if __name__ == "__main__":
    main()
