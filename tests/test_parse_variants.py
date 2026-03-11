"""
Tests for parse_variants.py
"""
import sys
from pathlib import Path

import pytest

# Allow importing from scripts/
sys.path.insert(0, str(Path(__file__).parent.parent / "scripts"))

from parse_variants import parse_variant, load_variants


class TestParseVariant:
    def test_valid_variant(self):
        v = parse_variant("R80C")
        assert v == {"wt_aa": "R", "position": 80, "mut_aa": "C", "raw": "R80C"}

    def test_lowercase_accepted(self):
        v = parse_variant("r80c")
        assert v["wt_aa"] == "R"
        assert v["mut_aa"] == "C"
        assert v["raw"] == "R80C"

    def test_invalid_format_raises(self):
        with pytest.raises(ValueError, match="Invalid variant format"):
            parse_variant("invalid")

    def test_synonymous_raises(self):
        with pytest.raises(ValueError, match="identical"):
            parse_variant("R80R")

    def test_unknown_wt_aa_raises(self):
        with pytest.raises(ValueError, match="Unknown wild-type"):
            parse_variant("X80C")

    def test_unknown_mut_aa_raises(self):
        with pytest.raises(ValueError, match="Unknown mutant"):
            parse_variant("R80X")

    def test_variant_with_large_position(self):
        v = parse_variant("K1234A")
        assert v["position"] == 1234

    def test_all_standard_amino_acids(self):
        # Spot-check a few valid pairs
        for wt, mut in [("A", "G"), ("L", "V"), ("N", "S"), ("W", "R")]:
            v = parse_variant(f"{wt}10{mut}")
            assert v["wt_aa"] == wt
            assert v["mut_aa"] == mut


class TestLoadVariants:
    def test_load_from_file(self, tmp_path):
        vf = tmp_path / "variants.txt"
        vf.write_text("R80C\nA123V\n")
        result = load_variants(vf)
        assert len(result) == 2
        assert result[0]["raw"] == "R80C"
        assert result[1]["raw"] == "A123V"

    def test_skip_comments_and_blanks(self, tmp_path):
        vf = tmp_path / "variants.txt"
        vf.write_text("# comment\n\nR80C\n\n# another\nA123V\n")
        result = load_variants(vf)
        assert len(result) == 2

    def test_file_not_found_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            load_variants(tmp_path / "nonexistent.txt")

    def test_invalid_line_raises_with_line_number(self, tmp_path):
        vf = tmp_path / "variants.txt"
        vf.write_text("R80C\nbadline\n")
        with pytest.raises(ValueError, match="Line 2"):
            load_variants(vf)

    def test_real_sample_file(self):
        """Load the bundled sample variants file."""
        sample = (
            Path(__file__).parent.parent / "data" / "variants" / "variants.txt"
        )
        result = load_variants(sample)
        assert len(result) == 5
        raws = [v["raw"] for v in result]
        assert "R80C" in raws
        assert "A123V" in raws
