"""
run_foldx.py – Execute FoldX BuildModel (and optionally AnalyseComplex).

FoldX must be installed separately and its path provided via ``--foldx``
or the ``FOLDX_BINARY`` environment variable.

Usage
-----
    # ΔΔG folding only
    python scripts/run_foldx.py \\
        --pdb          data/structures/protein.pdb \\
        --individual-list data/results/foldx_inputs/individual_list.txt \\
        --output-dir   data/results/foldx_outputs \\
        --foldx        /path/to/foldx

    # ΔΔG folding + binding (complex)
    python scripts/run_foldx.py \\
        --pdb          data/structures/complex.pdb \\
        --individual-list data/results/foldx_inputs/individual_list.txt \\
        --output-dir   data/results/foldx_outputs \\
        --foldx        /path/to/foldx \\
        --complex
"""

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path


def _find_foldx(foldx_path: str | None = None) -> str:
    """Resolve the FoldX binary path.

    Checks (in order):
    1. ``foldx_path`` argument
    2. ``FOLDX_BINARY`` environment variable
    3. ``foldx`` on PATH

    Raises
    ------
    FileNotFoundError
        If FoldX cannot be found.
    """
    candidates = []
    if foldx_path:
        candidates.append(foldx_path)
    env_path = os.environ.get("FOLDX_BINARY")
    if env_path:
        candidates.append(env_path)
    candidates.append("foldx")  # hope it's on PATH

    for c in candidates:
        if shutil.which(c) or Path(c).is_file():
            return c

    raise FileNotFoundError(
        "FoldX binary not found. "
        "Set FOLDX_BINARY environment variable or pass --foldx <path>."
    )


def run_build_model(
    pdb_path: str | Path,
    individual_list_path: str | Path,
    output_dir: str | Path,
    foldx_binary: str | None = None,
    number_of_runs: int = 5,
    dry_run: bool = False,
) -> subprocess.CompletedProcess | None:
    """Run FoldX BuildModel for each mutation in the individual list.

    Parameters
    ----------
    pdb_path : str or Path
        Path to the (repaired) PDB file.
    individual_list_path : str or Path
        Path to ``individual_list.txt``.
    output_dir : str or Path
        Directory where FoldX will write its output files.
    foldx_binary : str, optional
        Path to the FoldX executable.
    number_of_runs : int
        Number of stochastic runs for energy averaging (default: 5).
    dry_run : bool
        If ``True``, print the command but do not execute it.

    Returns
    -------
    subprocess.CompletedProcess or None
        Result of the subprocess call, or ``None`` if ``dry_run`` is ``True``.
    """
    pdb_path = Path(pdb_path).resolve()
    individual_list_path = Path(individual_list_path).resolve()
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    if not dry_run:
        foldx_bin = _find_foldx(foldx_binary)
    else:
        foldx_bin = foldx_binary or "foldx"

    cmd = [
        foldx_bin,
        "--command=BuildModel",
        f"--pdb={pdb_path.name}",
        f"--pdb-dir={pdb_path.parent}",
        f"--mutant-file={individual_list_path}",
        f"--output-dir={output_dir}",
        f"--numberOfRuns={number_of_runs}",
        "--out-pdb=false",
    ]

    print("Running FoldX BuildModel:")
    print("  " + " ".join(cmd))

    if dry_run:
        print("  [dry-run: not executing]")
        return None

    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        sys.stderr.write(f"FoldX BuildModel failed (exit {result.returncode}):\n")
        sys.stderr.write(result.stderr)
        result.check_returncode()

    return result


def run_analyse_complex(
    pdb_path: str | Path,
    output_dir: str | Path,
    chain_of_interest: str = "A",
    foldx_binary: str | None = None,
    dry_run: bool = False,
) -> subprocess.CompletedProcess | None:
    """Run FoldX AnalyseComplex to compute ΔG binding.

    Should be called on both the wild-type and each mutant PDB to later
    compute ΔΔG binding = ΔG(mutant) – ΔG(WT).

    Parameters
    ----------
    pdb_path : str or Path
        PDB file to analyse.
    output_dir : str or Path
        Directory for output files.
    chain_of_interest : str
        Chain of the protein of interest (default: ``"A"``).
    foldx_binary : str, optional
    dry_run : bool

    Returns
    -------
    subprocess.CompletedProcess or None
    """
    pdb_path = Path(pdb_path).resolve()
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    if not dry_run:
        foldx_bin = _find_foldx(foldx_binary)
    else:
        foldx_bin = foldx_binary or "foldx"

    cmd = [
        foldx_bin,
        "--command=AnalyseComplex",
        f"--pdb={pdb_path.name}",
        f"--pdb-dir={pdb_path.parent}",
        f"--analyseComplexChains={chain_of_interest}",
        f"--output-dir={output_dir}",
    ]

    print("Running FoldX AnalyseComplex:")
    print("  " + " ".join(cmd))

    if dry_run:
        print("  [dry-run: not executing]")
        return None

    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        sys.stderr.write(f"FoldX AnalyseComplex failed (exit {result.returncode}):\n")
        sys.stderr.write(result.stderr)
        result.check_returncode()

    return result


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Run FoldX BuildModel (and optionally AnalyseComplex) for each mutation."
    )
    parser.add_argument("--pdb", required=True, metavar="PDB",
                        help="Path to the PDB structure file.")
    parser.add_argument("--individual-list", required=True, metavar="FILE",
                        help="Path to FoldX individual_list.txt.")
    parser.add_argument("--output-dir", default="data/results/foldx_outputs",
                        metavar="DIR", help="Directory for FoldX output files.")
    parser.add_argument("--foldx", default=None, metavar="BINARY",
                        help="Path to the FoldX executable.")
    parser.add_argument("--runs", type=int, default=5, metavar="N",
                        help="Number of stochastic FoldX runs (default: 5).")
    parser.add_argument("--complex", action="store_true",
                        help="Also run AnalyseComplex to compute ΔΔG binding.")
    parser.add_argument("--chain", default="A", metavar="CHAIN",
                        help="Chain of interest for AnalyseComplex (default: A).")
    parser.add_argument("--dry-run", action="store_true",
                        help="Print FoldX commands without executing them.")
    args = parser.parse_args(argv)

    run_build_model(
        pdb_path=args.pdb,
        individual_list_path=args.individual_list,
        output_dir=args.output_dir,
        foldx_binary=args.foldx,
        number_of_runs=args.runs,
        dry_run=args.dry_run,
    )

    if args.complex:
        run_analyse_complex(
            pdb_path=args.pdb,
            output_dir=args.output_dir,
            chain_of_interest=args.chain,
            foldx_binary=args.foldx,
            dry_run=args.dry_run,
        )


if __name__ == "__main__":
    main()
