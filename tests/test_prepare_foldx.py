"""
Tests for prepare_foldx.py
"""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent / "scripts"))

from prepare_foldx import format_foldx_mutation, write_individual_list, prepare_foldx_inputs


def _variant(wt, pos, mut, chain="A"):
    return {"wt_aa": wt, "position": pos, "mut_aa": mut, "raw": f"{wt}{pos}{mut}", "chain": chain}


class TestFormatFoldxMutation:
    def test_basic(self):
        v = _variant("R", 80, "C")
        assert format_foldx_mutation(v) == "AR80C"

    def test_chain_b(self):
        v = _variant("A", 10, "G", chain="B")
        assert format_foldx_mutation(v) == "BA10G"

    def test_default_chain_a(self):
        v = {"wt_aa": "L", "position": 200, "mut_aa": "P"}
        assert format_foldx_mutation(v) == "AL200P"


class TestWriteIndividualList:
    def test_one_per_line(self, tmp_path):
        variants = [_variant("R", 80, "C"), _variant("A", 123, "V")]
        out = tmp_path / "individual_list.txt"
        write_individual_list(variants, out, one_per_line=True)
        lines = out.read_text().splitlines()
        assert lines == ["AR80C;", "AA123V;"]

    def test_combined_single_line(self, tmp_path):
        variants = [_variant("R", 80, "C"), _variant("A", 123, "V")]
        out = tmp_path / "individual_list.txt"
        write_individual_list(variants, out, one_per_line=False)
        content = out.read_text().strip()
        assert content == "AR80C,AA123V;"

    def test_creates_parent_dirs(self, tmp_path):
        out = tmp_path / "nested" / "dir" / "individual_list.txt"
        write_individual_list([_variant("R", 80, "C")], out)
        assert out.exists()

    def test_single_variant(self, tmp_path):
        out = tmp_path / "individual_list.txt"
        write_individual_list([_variant("G", 45, "D")], out)
        assert out.read_text().strip() == "AG45D;"


class TestPrepareFoldxInputs:
    def test_creates_individual_list(self, tmp_path):
        variants = [_variant("R", 80, "C"), _variant("L", 200, "P")]
        result = prepare_foldx_inputs(variants, tmp_path / "foldx_in")
        assert result["individual_list"].exists()
        content = result["individual_list"].read_text()
        assert "AR80C;" in content
        assert "AL200P;" in content
