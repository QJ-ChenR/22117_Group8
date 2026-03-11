"""
prepare_foldx.py – Generate FoldX input files for BuildModel runs.

FoldX expects an ``individual_list.txt`` file where each line lists the
mutations for one run, in the format::

    <chain><WT_AA><position><MUT_AA>;

For example, mutation R80C on chain A is written as ``AR80C;``.
Multiple mutations on the same line are separated by commas, e.g.::

    AR80C,AD45G;

Usage
-----
    python scripts/prepare_foldx.py \\
        --variants   data/variants/variants.txt \\
        --chain      A \\
        --output-dir data/results/foldx_inputs
"""

import argparse
from pathlib import Path

from parse_variants import load_variants


def format_foldx_mutation(variant: dict) -> str:
    """Format a single variant as a FoldX mutation string.

    Parameters
    ----------
    variant : dict
        Variant dict (from :func:`parse_variants.load_variants` or the mapped
        dicts produced by :mod:`map_to_structure`).  Must contain keys
        ``chain``, ``wt_aa``, ``position``, ``mut_aa``.  If ``chain`` is
        absent it defaults to ``"A"``.

    Returns
    -------
    str
        FoldX mutation token, e.g. ``"AR80C"``.
    """
    chain = variant.get("chain", "A")
    return f"{chain}{variant['wt_aa']}{variant['position']}{variant['mut_aa']}"


def write_individual_list(
    variants: list[dict],
    output_path: str | Path,
    one_per_line: bool = True,
) -> Path:
    """Write a FoldX ``individual_list.txt`` file.

    Parameters
    ----------
    variants : list of dict
        Variants to include.  Each dict must have ``wt_aa``, ``position``,
        ``mut_aa`` and optionally ``chain`` (default ``"A"``).
    output_path : str or Path
        Destination file path (will be created, overwriting if present).
    one_per_line : bool
        If ``True`` (default) each variant is written on its own line so that
        FoldX evaluates them independently.  If ``False`` all mutations are
        written on a single line (simultaneous multi-point mutagenesis).

    Returns
    -------
    Path
        Absolute path to the written file.
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with output_path.open("w") as fh:
        if one_per_line:
            for v in variants:
                fh.write(format_foldx_mutation(v) + ";\n")
        else:
            mutations = ",".join(format_foldx_mutation(v) for v in variants)
            fh.write(mutations + ";\n")

    return output_path.resolve()


def prepare_foldx_inputs(
    variants: list[dict],
    output_dir: str | Path,
) -> dict:
    """Prepare the directory and files needed for a FoldX BuildModel run.

    Creates ``<output_dir>/individual_list.txt``.

    Parameters
    ----------
    variants : list of dict
    output_dir : str or Path

    Returns
    -------
    dict
        ``{"individual_list": Path, "variants": list}``
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    individual_list_path = write_individual_list(
        variants, output_dir / "individual_list.txt"
    )

    return {
        "individual_list": individual_list_path,
        "variants": variants,
    }


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Prepare FoldX individual_list.txt from a variants file."
    )
    parser.add_argument("--variants", required=True, metavar="FILE",
                        help="Path to the variants file.")
    parser.add_argument("--chain", default="A", metavar="CHAIN",
                        help="Chain ID to assign to all mutations (default: A).")
    parser.add_argument("--output-dir", default="data/results/foldx_inputs",
                        metavar="DIR",
                        help="Directory where individual_list.txt will be written.")
    parser.add_argument("--combined", action="store_true",
                        help="Write all mutations on a single line (simultaneous mutagenesis).")
    args = parser.parse_args(argv)

    variants = load_variants(args.variants)
    # Inject chain into each variant dict
    for v in variants:
        v.setdefault("chain", args.chain)

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = write_individual_list(
        variants, out_dir / "individual_list.txt", one_per_line=not args.combined
    )

    print(f"Wrote FoldX individual_list.txt → {out_path}")
    with out_path.open() as fh:
        print("\nContents:")
        print(fh.read())

    return out_path


if __name__ == "__main__":
    main()
