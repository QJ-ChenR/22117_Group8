"""
summarize_results.py – Produce a summary CSV table of ΔΔG values.

Merges ΔΔG folding and (optionally) ΔΔG binding results into a single
pandas DataFrame and writes a CSV file.

Usage
-----
    python scripts/summarize_results.py \\
        --output-dir  data/results/foldx_outputs \\
        --pdb-name    protein \\
        --summary-out data/results/summary.csv

    # Including binding (complex):
    python scripts/summarize_results.py \\
        --output-dir  data/results/foldx_outputs \\
        --pdb-name    complex \\
        --summary-out data/results/summary.csv \\
        --complex
"""

import argparse
from pathlib import Path

try:
    import pandas as pd
except ImportError as exc:  # pragma: no cover
    raise ImportError(
        "pandas is required. Install it with: pip install pandas"
    ) from exc

from parse_foldx_output import load_foldx_results


def build_summary_table(
    ddg_folding: list[dict],
    ddg_binding: list[dict] | None = None,
) -> "pd.DataFrame":
    """Build a summary DataFrame.

    Parameters
    ----------
    ddg_folding : list of dict
        Each dict: ``{"mutation": str, "ddG_folding": float}``
    ddg_binding : list of dict, optional
        Each dict: ``{"mutation": str, "ddG_binding": float}``

    Returns
    -------
    pd.DataFrame
        Columns: ``mutation``, ``ddG_folding``, (``ddG_binding`` if provided).
        Rows are sorted by ``mutation``.
    """
    df = pd.DataFrame(ddg_folding)
    if df.empty:
        df = pd.DataFrame(columns=["mutation", "ddG_folding"])

    if ddg_binding:
        df_binding = pd.DataFrame(ddg_binding)
        df = df.merge(df_binding, on="mutation", how="left")

    df = df.sort_values("mutation").reset_index(drop=True)
    return df


def summarize(
    output_dir: str | Path,
    pdb_name: str,
    summary_out: str | Path,
    is_complex: bool = False,
) -> "pd.DataFrame":
    """Load FoldX results, build table, and write CSV.

    Parameters
    ----------
    output_dir : str or Path
    pdb_name : str
    summary_out : str or Path
    is_complex : bool

    Returns
    -------
    pd.DataFrame
    """
    results = load_foldx_results(output_dir, pdb_name, is_complex)

    df = build_summary_table(
        results["ddg_folding"],
        results["ddg_binding"] if is_complex else None,
    )

    summary_out = Path(summary_out)
    summary_out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(summary_out, index=False)
    print(f"Summary written → {summary_out}")

    return df


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Summarise FoldX ΔΔG results into a CSV table."
    )
    parser.add_argument("--output-dir", required=True, metavar="DIR",
                        help="Directory containing FoldX output files.")
    parser.add_argument("--pdb-name", required=True, metavar="NAME",
                        help="PDB base name without extension (e.g. 'protein').")
    parser.add_argument("--summary-out", default="data/results/summary.csv",
                        metavar="FILE",
                        help="Path for the output CSV file (default: data/results/summary.csv).")
    parser.add_argument("--complex", action="store_true",
                        help="Also include ΔΔG binding from AnalyseComplex.")
    args = parser.parse_args(argv)

    df = summarize(args.output_dir, args.pdb_name, args.summary_out, args.complex)

    print("\nSummary table:")
    print(df.to_string(index=False))

    return df


if __name__ == "__main__":
    main()
