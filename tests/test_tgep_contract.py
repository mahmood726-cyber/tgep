import shutil
import subprocess
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]


def test_tgep_r_exports_core_functions_and_guards():
    text = (ROOT / "R" / "TGEP.R").read_text(encoding="utf-8")
    required = [
        "grma_guard_core <- function",
        "wrd_guard_core <- function",
        "swa_guard_core <- function",
        "tgep_meta <- function",
        "compare_tgep <- function",
        "plot.tgep <- function",
    ]
    missing = [marker for marker in required if marker not in text]
    assert missing == []


def test_existing_r_suite_documents_edge_cases():
    text = (ROOT / "test_tgep.R").read_text(encoding="utf-8")
    required = [
        "K=2 fallback",
        "K=1 fallback",
        "All-identical studies",
        "Input validation",
        "Publication bias scenario",
    ]
    missing = [marker for marker in required if marker not in text]
    assert missing == []


def test_r_suite_passes_when_rscript_is_available():
    rscript = shutil.which("Rscript")
    if rscript is None:
        pytest.skip("Rscript is not installed in this environment")
    subprocess.run([rscript, "test_tgep.R"], cwd=ROOT, check=True)
