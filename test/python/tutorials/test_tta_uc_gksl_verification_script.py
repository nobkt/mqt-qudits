"""Tests for TTA-UC GKSL verification runner script."""

from __future__ import annotations

import json
import sys
from pathlib import Path

# Add tutorials to path
TUTORIALS_PATH = Path(__file__).parent.parent.parent.parent / "tutorials"
sys.path.insert(0, str(TUTORIALS_PATH))

from run_tta_uc_gksl_verification import run_verification_workflow


def test_run_verification_workflow_writes_reports(tmp_path: Path) -> None:
    json_path, markdown_path, report = run_verification_workflow(
        output_dir=tmp_path,
        scenarios=("classical", "qudit"),
        t_max=2.0,
        n_steps=2,
        initial_state="edge_triplet",
    )

    assert json_path.exists()
    assert markdown_path.exists()
    assert report["success"]
    assert report["errors"] == []
    assert [entry["scenario"] for entry in report["results"]] == ["classical", "qudit"]

    loaded = json.loads(json_path.read_text(encoding="utf-8"))
    assert loaded["success"]
    assert len(loaded["results"]) == 2
