# 22117_Group8 – MAVISp-style ΔΔG Pipeline

A modular, command-line Python pipeline that reproduces a simplified
[MAVISp](https://github.com/ELELAB/MAVISp)-style workflow for analysing the
effect of missense mutations on a human protein structure.

---

## Overview

Given a list of missense variants (e.g. `R80C`) and a PDB/AlphaFold structure,
the pipeline:

1. Parses and validates each missense variant.
2. Maps variants onto the PDB structure (checks wild-type residue identity).
3. Prepares FoldX `individual_list.txt` input files.
4. Runs **FoldX BuildModel** to compute **ΔΔG folding**.
5. Optionally runs **FoldX AnalyseComplex** to compute **ΔΔG binding** (for
   protein complexes).
6. Parses FoldX output files.
7. Produces a summary CSV table with columns `mutation`, `ddG_folding`
   (and `ddG_binding` if applicable).

---

## Repository Layout

```
data/
  structures/        PDB or AlphaFold structure files
  variants/          Missense variant lists (one variant per line)
  results/           Pipeline outputs (FoldX inputs/outputs, summary CSV)
scripts/
  parse_variants.py      Step 1 – parse variant strings
  map_to_structure.py    Step 2 – validate variants against PDB
  prepare_foldx.py       Step 3 – write FoldX individual_list.txt
  run_foldx.py           Step 4 – execute FoldX
  parse_foldx_output.py  Step 5 – parse FoldX output files
  summarize_results.py   Step 6 – build summary CSV
  pipeline.py            Main orchestrator (runs all steps)
tests/                 pytest test suite
requirements.txt
```

---

## Quick Start

### 1. Install dependencies

```bash
pip install -r requirements.txt
```

> **Note:** FoldX must be installed separately.
> Download from [https://foldxsuite.crg.eu/](https://foldxsuite.crg.eu/) and
> either add it to your PATH or set the `FOLDX_BINARY` environment variable.

### 2. Prepare input files

* **Variants file** (`data/variants/variants.txt`):

  ```
  # One missense variant per line (WT_AA + position + MUT_AA)
  R80C
  A123V
  G45D
  ```

* **Structure file** (`data/structures/protein.pdb`):
  A PDB or AlphaFold `.pdb` file for your protein of interest.

### 3. Run the pipeline

```bash
# Minimal run (ΔΔG folding only)
python scripts/pipeline.py \
    --variants  data/variants/variants.txt \
    --structure data/structures/protein.pdb \
    --chain     A \
    --output-dir data/results \
    --foldx     /path/to/foldx

# Complex run (also compute ΔΔG binding)
python scripts/pipeline.py \
    --variants  data/variants/variants.txt \
    --structure data/structures/complex.pdb \
    --chain     A \
    --output-dir data/results \
    --foldx     /path/to/foldx \
    --complex

# Dry-run (validate inputs and print FoldX commands without executing)
python scripts/pipeline.py \
    --variants  data/variants/variants.txt \
    --structure data/structures/protein.pdb \
    --dry-run
```

The final summary CSV (`data/results/summary.csv`) contains:

| mutation | ddG_folding | ddG_binding |
|----------|-------------|-------------|
| AR80C    | 2.50        | 1.20        |
| AA123V   | -1.00       | —           |

---

## Running Individual Steps

Each script is independently runnable:

```bash
# Parse variants
python scripts/parse_variants.py --variants data/variants/variants.txt

# Map to structure
python scripts/map_to_structure.py \
    --structure data/structures/protein.pdb \
    --variants  data/variants/variants.txt \
    --chain     A

# Prepare FoldX inputs
python scripts/prepare_foldx.py \
    --variants   data/variants/variants.txt \
    --chain      A \
    --output-dir data/results/foldx_inputs

# Run FoldX
python scripts/run_foldx.py \
    --pdb              data/structures/protein.pdb \
    --individual-list  data/results/foldx_inputs/individual_list.txt \
    --output-dir       data/results/foldx_outputs \
    --foldx            /path/to/foldx

# Parse FoldX outputs
python scripts/parse_foldx_output.py \
    --output-dir data/results/foldx_outputs \
    --pdb-name   protein

# Build summary CSV
python scripts/summarize_results.py \
    --output-dir  data/results/foldx_outputs \
    --pdb-name    protein \
    --summary-out data/results/summary.csv
```

---

## Running Tests

```bash
pip install pytest
pytest tests/ -v
```

---

## FoldX Notes

* The pipeline uses **FoldX BuildModel** (ΔΔG folding) and optionally
  **AnalyseComplex** (ΔΔG binding).
* FoldX expects a "repaired" PDB. For best results, run
  `FoldX --command=RepairPDB` on your structure before using this pipeline.
* `numberOfRuns` (default: 5) controls how many stochastic runs are averaged.

---

## Dependencies

| Package    | Version   | Purpose                      |
|------------|-----------|------------------------------|
| biopython  | ≥ 1.81    | PDB parsing                  |
| pandas     | ≥ 2.0     | Data manipulation / CSV output|
| pytest     | any       | Running tests                |
