"""
Tests for parse_foldx_output.py
"""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent / "scripts"))

from parse_foldx_output import (
    parse_dif_buildmodel,
    parse_average_buildmodel,
    extract_ddg_folding,
    parse_summary_ac,
    compute_ddg_binding,
    load_foldx_results,
    _parse_tsv_block,
)


# ── Helpers ───────────────────────────────────────────────────────────────────

def _write_dif_buildmodel(tmp_path: Path, rows: list[str]) -> Path:
    """Write a mock Dif_BuildModel_protein.fxout file."""
    header = "Pdb\ttotal energy\tBackbone Hbond\tSidechain Hbond\tVan der Waals\n"
    content = header + "\n".join(rows) + "\n"
    p = tmp_path / "Dif_BuildModel_protein.fxout"
    p.write_text(content)
    return p


def _write_summary_ac(tmp_path: Path, rows: list[str], filename: str) -> Path:
    header = "Pdb\tGroup 1\tGroup 2\tInteraction Energy\tBackbone Hbond\n"
    content = header + "\n".join(rows) + "\n"
    p = tmp_path / filename
    p.write_text(content)
    return p


# ── Tests: _parse_tsv_block ───────────────────────────────────────────────────

class TestParseTsvBlock:
    def test_basic(self):
        lines = ["Col1\tCol2\n", "a\t1\n", "b\t2\n"]
        result = _parse_tsv_block(lines)
        assert result == [{"Col1": "a", "Col2": "1"}, {"Col1": "b", "Col2": "2"}]

    def test_skips_comment_lines(self):
        lines = ["Col1\tCol2\n", "// comment\n", "a\t1\n"]
        result = _parse_tsv_block(lines)
        assert len(result) == 1

    def test_empty_returns_empty(self):
        assert _parse_tsv_block([]) == []


# ── Tests: parse_dif_buildmodel ───────────────────────────────────────────────

class TestParseDifBuildmodel:
    def test_parses_rows(self, tmp_path):
        f = _write_dif_buildmodel(tmp_path, [
            "AR80C_1.pdb\t2.50\t0.1\t0.2\t0.3",
            "AA123V_1.pdb\t-1.00\t0.0\t0.0\t0.0",
        ])
        rows = parse_dif_buildmodel(f)
        assert len(rows) == 2
        assert rows[0]["Pdb"] == "AR80C_1.pdb"
        assert rows[0]["total energy"] == "2.50"

    def test_file_not_found(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            parse_dif_buildmodel(tmp_path / "missing.fxout")


# ── Tests: extract_ddg_folding ────────────────────────────────────────────────

class TestExtractDdgFolding:
    def test_extracts_values(self, tmp_path):
        f = _write_dif_buildmodel(tmp_path, [
            "AR80C_1.pdb\t2.50\t0\t0\t0",
            "AA123V_1.pdb\t-1.00\t0\t0\t0",
        ])
        rows = parse_dif_buildmodel(f)
        result = extract_ddg_folding(rows)
        assert result[0]["mutation"] == "AR80C"
        assert result[0]["ddG_folding"] == pytest.approx(2.50)
        assert result[1]["ddG_folding"] == pytest.approx(-1.00)

    def test_nan_for_bad_value(self, tmp_path):
        f = _write_dif_buildmodel(tmp_path, ["AR80C_1.pdb\tN/A\t0\t0\t0"])
        rows = parse_dif_buildmodel(f)
        result = extract_ddg_folding(rows)
        import math
        assert math.isnan(result[0]["ddG_folding"])


# ── Tests: compute_ddg_binding ────────────────────────────────────────────────

class TestComputeDdgBinding:
    def test_basic_subtraction(self, tmp_path):
        wt = [{"Pdb": "wt.pdb", "Interaction Energy": "-10.0"}]
        mut = [{"Pdb": "AR80C_1.pdb", "Interaction Energy": "-8.0"}]
        result = compute_ddg_binding(wt, mut)
        assert result[0]["ddG_binding"] == pytest.approx(2.0)

    def test_empty_wt_returns_empty(self):
        result = compute_ddg_binding([], [{"Pdb": "AR80C_1.pdb", "Interaction Energy": "1.0"}])
        assert result == []

    def test_stabilising_mutation(self, tmp_path):
        wt = [{"Pdb": "wt.pdb", "Interaction Energy": "-10.0"}]
        mut = [{"Pdb": "AA123V_1.pdb", "Interaction Energy": "-12.0"}]
        result = compute_ddg_binding(wt, mut)
        assert result[0]["ddG_binding"] == pytest.approx(-2.0)


# ── Tests: load_foldx_results ─────────────────────────────────────────────────

class TestLoadFoldxResults:
    def test_loads_folding_only(self, tmp_path):
        _write_dif_buildmodel(tmp_path, [
            "AR80C_1.pdb\t1.5\t0\t0\t0",
        ])
        result = load_foldx_results(tmp_path, "protein", is_complex=False)
        assert len(result["ddg_folding"]) == 1
        assert result["ddg_binding"] == []

    def test_no_dif_file_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            load_foldx_results(tmp_path, "nonexistent")
