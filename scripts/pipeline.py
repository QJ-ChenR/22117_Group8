"""
pipeline.py – Main CLI orchestrator for the MAVISp-style ΔΔG pipeline.

Ties together all steps:
  1. Parse missense variants
  2. Load and validate against PDB structure
  3. Prepare FoldX inputs
  4. Run FoldX BuildModel (and optionally AnalyseComplex)
  5. Parse FoldX outputs
  6. Produce summary CSV

Usage
-----
    python scripts/pipeline.py \\
        --variants   data/variants/variants.txt \\
        --structure  data/structures/protein.pdb \\
        --chain      A \\
        --output-dir data/results \\
        --foldx      /path/to/foldx

    # Full complex run (also compute ΔΔG binding):
    python scripts/pipeline.py \\
        --variants   data/variants/variants.txt \\
        --structure  data/structures/complex.pdb \\
        --chain      A \\
        --output-dir data/results \\
        --foldx      /path/to/foldx \\
        --complex

    # Dry-run (skip FoldX execution):
    python scripts/pipeline.py \\
        --variants   data/variants/variants.txt \\
        --structure  data/structures/protein.pdb \\
        --output-dir data/results \\
        --dry-run
"""

import argparse
import sys
from pathlib import Path


def _add_scripts_to_path():
    """Ensure the scripts/ directory is on sys.path when invoked directly."""
    scripts_dir = Path(__file__).parent.resolve()
    if str(scripts_dir) not in sys.path:
        sys.path.insert(0, str(scripts_dir))


_add_scripts_to_path()

from parse_variants import load_variants  # noqa: E402
from map_to_structure import load_structure, map_variants_to_structure  # noqa: E402
from prepare_foldx import prepare_foldx_inputs  # noqa: E402
from run_foldx import run_build_model, run_analyse_complex  # noqa: E402
from parse_foldx_output import load_foldx_results  # noqa: E402
from summarize_results import build_summary_table  # noqa: E402


def run_pipeline(
    variants_file: str | Path,
    structure_file: str | Path,
    output_dir: str | Path,
    chain: str = "A",
    foldx_binary: str | None = None,
    number_of_runs: int = 5,
    is_complex: bool = False,
    dry_run: bool = False,
    skip_foldx: bool = False,
) -> dict:
    """Execute the full ΔΔG pipeline.

    Parameters
    ----------
    variants_file : str or Path
        Path to the variants file (one variant per line, e.g. "R80C").
    structure_file : str or Path
        Path to the PDB structure file.
    output_dir : str or Path
        Root directory for pipeline outputs.
    chain : str
        Chain ID to use for mutation mapping and FoldX (default: "A").
    foldx_binary : str, optional
        Path to the FoldX executable.
    number_of_runs : int
        Number of stochastic FoldX BuildModel runs (default: 5).
    is_complex : bool
        If ``True``, also run AnalyseComplex to compute ΔΔG binding.
    dry_run : bool
        Print FoldX commands without executing them.
    skip_foldx : bool
        Skip FoldX steps entirely (useful for testing the upstream steps).

    Returns
    -------
    dict
        ``{"variants": list, "mapped": list, "summary": pd.DataFrame}``
    """
    output_dir = Path(output_dir)
    foldx_inputs_dir = output_dir / "foldx_inputs"
    foldx_outputs_dir = output_dir / "foldx_outputs"
    summary_path = output_dir / "summary.csv"
    pdb_name = Path(structure_file).stem

    # ── Step 1: Parse variants ────────────────────────────────────────────────
    print("=" * 60)
    print("Step 1: Parsing variants")
    print("=" * 60)
    variants = load_variants(variants_file)
    print(f"  Loaded {len(variants)} variant(s).")
    for v in variants:
        print(f"    {v['raw']}")

    # ── Step 2: Map to structure ──────────────────────────────────────────────
    print("\n" + "=" * 60)
    print("Step 2: Mapping variants to structure")
    print("=" * 60)
    structure = load_structure(structure_file)
    mapped = map_variants_to_structure(structure, variants, chain)

    ok = [v for v in mapped if v["mapped"]]
    failed = [v for v in mapped if not v["mapped"]]
    print(f"  {len(ok)}/{len(mapped)} variant(s) mapped successfully.")
    for v in failed:
        print(f"  WARNING: {v['raw']} – {v['warning']}")

    if not ok:
        print("  No variants could be mapped. Aborting.")
        return {"variants": variants, "mapped": mapped, "summary": None}

    # Use only successfully mapped variants downstream
    variants_to_use = ok

    # ── Step 3: Prepare FoldX inputs ──────────────────────────────────────────
    print("\n" + "=" * 60)
    print("Step 3: Preparing FoldX inputs")
    print("=" * 60)
    foldx_prep = prepare_foldx_inputs(variants_to_use, foldx_inputs_dir)
    print(f"  individual_list.txt → {foldx_prep['individual_list']}")

    if skip_foldx or dry_run:
        action = "Skipping" if skip_foldx else "Dry-run"
        print(f"\n  {action}: FoldX steps will not be executed.")
        return {
            "variants": variants,
            "mapped": mapped,
            "summary": None,
            "individual_list": foldx_prep["individual_list"],
        }

    # ── Step 4: Run FoldX ─────────────────────────────────────────────────────
    print("\n" + "=" * 60)
    print("Step 4: Running FoldX BuildModel")
    print("=" * 60)
    run_build_model(
        pdb_path=structure_file,
        individual_list_path=foldx_prep["individual_list"],
        output_dir=foldx_outputs_dir,
        foldx_binary=foldx_binary,
        number_of_runs=number_of_runs,
        dry_run=dry_run,
    )

    if is_complex:
        print("\n" + "=" * 60)
        print("Step 4b: Running FoldX AnalyseComplex")
        print("=" * 60)
        run_analyse_complex(
            pdb_path=structure_file,
            output_dir=foldx_outputs_dir,
            chain_of_interest=chain,
            foldx_binary=foldx_binary,
            dry_run=dry_run,
        )

    # ── Step 5: Parse FoldX outputs ───────────────────────────────────────────
    print("\n" + "=" * 60)
    print("Step 5: Parsing FoldX outputs")
    print("=" * 60)
    foldx_results = load_foldx_results(foldx_outputs_dir, pdb_name, is_complex)

    # ── Step 6: Summarise ─────────────────────────────────────────────────────
    print("\n" + "=" * 60)
    print("Step 6: Building summary table")
    print("=" * 60)
    df = build_summary_table(
        foldx_results["ddg_folding"],
        foldx_results["ddg_binding"] if is_complex else None,
    )
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(summary_path, index=False)
    print(f"  Summary CSV → {summary_path}")
    print("\n" + df.to_string(index=False))

    return {"variants": variants, "mapped": mapped, "summary": df}


def main(argv=None):
    parser = argparse.ArgumentParser(
        description=(
            "MAVISp-style pipeline: compute ΔΔG folding (and optionally "
            "ΔΔG binding) for a list of missense variants."
        )
    )
    parser.add_argument("--variants", required=True, metavar="FILE",
                        help="Path to the variants file (one variant per line).")
    parser.add_argument("--structure", required=True, metavar="PDB",
                        help="Path to the PDB structure file.")
    parser.add_argument("--chain", default="A", metavar="CHAIN",
                        help="Chain ID to use (default: A).")
    parser.add_argument("--output-dir", default="data/results", metavar="DIR",
                        help="Root directory for all output files (default: data/results).")
    parser.add_argument("--foldx", default=None, metavar="BINARY",
                        help="Path to the FoldX executable (or set FOLDX_BINARY env var).")
    parser.add_argument("--runs", type=int, default=5, metavar="N",
                        help="Number of stochastic FoldX BuildModel runs (default: 5).")
    parser.add_argument("--complex", action="store_true",
                        help="Run AnalyseComplex to compute ΔΔG binding for a protein complex.")
    parser.add_argument("--dry-run", action="store_true",
                        help="Print FoldX commands without executing them.")
    parser.add_argument("--skip-foldx", action="store_true",
                        help="Skip FoldX execution entirely (parse + prepare steps only).")
    args = parser.parse_args(argv)

    run_pipeline(
        variants_file=args.variants,
        structure_file=args.structure,
        output_dir=args.output_dir,
        chain=args.chain,
        foldx_binary=args.foldx,
        number_of_runs=args.runs,
        is_complex=args.complex,
        dry_run=args.dry_run,
        skip_foldx=args.skip_foldx,
    )


if __name__ == "__main__":
    main()
