"""
Tests for map_to_structure.py
"""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent / "scripts"))

from map_to_structure import load_structure, get_residue_aa, map_variants_to_structure
from parse_variants import parse_variant

# Path to the bundled example PDB
EXAMPLE_PDB = Path(__file__).parent.parent / "data" / "structures" / "example_protein.pdb"


@pytest.fixture
def structure():
    return load_structure(EXAMPLE_PDB)


class TestLoadStructure:
    def test_loads_without_error(self, structure):
        assert structure is not None

    def test_raises_for_missing_file(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            load_structure(tmp_path / "missing.pdb")


class TestGetResidueAA:
    def test_known_residue(self, structure):
        aa = get_residue_aa(structure, "A", 80)
        assert aa == "R"  # ARG at position 80

    def test_known_residue_gly(self, structure):
        aa = get_residue_aa(structure, "A", 45)
        assert aa == "G"

    def test_missing_position_returns_none(self, structure):
        aa = get_residue_aa(structure, "A", 9999)
        assert aa is None

    def test_wrong_chain_returns_none(self, structure):
        aa = get_residue_aa(structure, "B", 80)
        assert aa is None


class TestMapVariantsToStructure:
    def _v(self, raw, chain="A"):
        v = parse_variant(raw)
        v["chain"] = chain
        return v

    def test_successful_mapping(self, structure):
        variants = [self._v("R80C"), self._v("A123V"), self._v("G45D")]
        mapped = map_variants_to_structure(structure, variants, "A")
        assert all(m["mapped"] for m in mapped)
        assert all(m["warning"] is None for m in mapped)

    def test_wt_mismatch_not_mapped(self, structure):
        # Position 80 is ARG (R), providing V as WT → mismatch
        variants = [parse_variant("V80C")]
        mapped = map_variants_to_structure(structure, variants, "A")
        assert not mapped[0]["mapped"]
        assert "mismatch" in mapped[0]["warning"]

    def test_missing_residue_not_mapped(self, structure):
        variants = [parse_variant("R9999C")]
        mapped = map_variants_to_structure(structure, variants, "A")
        assert not mapped[0]["mapped"]
        assert "not found" in mapped[0]["warning"]

    def test_mixed_results(self, structure):
        variants = [parse_variant("R80C"), parse_variant("R9999C")]
        mapped = map_variants_to_structure(structure, variants, "A")
        assert mapped[0]["mapped"] is True
        assert mapped[1]["mapped"] is False
