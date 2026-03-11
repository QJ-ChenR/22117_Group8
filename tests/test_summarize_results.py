"""
Tests for summarize_results.py
"""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent / "scripts"))

from summarize_results import build_summary_table


class TestBuildSummaryTable:
    def test_folding_only(self):
        ddg = [
            {"mutation": "AR80C", "ddG_folding": 2.5},
            {"mutation": "AA123V", "ddG_folding": -1.0},
        ]
        df = build_summary_table(ddg)
        assert list(df.columns) == ["mutation", "ddG_folding"]
        assert len(df) == 2

    def test_sorted_by_mutation(self):
        ddg = [
            {"mutation": "AR80C", "ddG_folding": 2.5},
            {"mutation": "AA10G", "ddG_folding": 0.1},
        ]
        df = build_summary_table(ddg)
        assert df.iloc[0]["mutation"] == "AA10G"
        assert df.iloc[1]["mutation"] == "AR80C"

    def test_with_binding(self):
        ddg_folding = [
            {"mutation": "AR80C", "ddG_folding": 2.5},
        ]
        ddg_binding = [
            {"mutation": "AR80C", "ddG_binding": 1.2},
        ]
        df = build_summary_table(ddg_folding, ddg_binding)
        assert "ddG_binding" in df.columns
        assert df.iloc[0]["ddG_binding"] == pytest.approx(1.2)

    def test_empty_input(self):
        df = build_summary_table([])
        assert len(df) == 0
        assert "mutation" in df.columns

    def test_binding_merge_preserves_all_mutations(self):
        ddg_folding = [
            {"mutation": "AR80C", "ddG_folding": 2.5},
            {"mutation": "AA123V", "ddG_folding": 1.0},
        ]
        # Only one of the two mutations has binding data
        ddg_binding = [{"mutation": "AR80C", "ddG_binding": 0.5}]
        df = build_summary_table(ddg_folding, ddg_binding)
        assert len(df) == 2
        aa123v_row = df[df["mutation"] == "AA123V"]
        import math
        assert math.isnan(aa123v_row.iloc[0]["ddG_binding"])

    def test_write_csv(self, tmp_path):
        ddg = [{"mutation": "AR80C", "ddG_folding": 2.5}]
        df = build_summary_table(ddg)
        csv_path = tmp_path / "summary.csv"
        df.to_csv(csv_path, index=False)
        import pandas as pd
        loaded = pd.read_csv(csv_path)
        assert loaded.iloc[0]["mutation"] == "AR80C"
        assert loaded.iloc[0]["ddG_folding"] == pytest.approx(2.5)
