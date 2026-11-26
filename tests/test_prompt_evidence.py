"""prompt-evidence: masked OpenGenes / SynergyAge prompt parsers (no DB)."""
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1] / "skills" / "longevity-analyst"
sys.path.insert(0, str(ROOT / "scripts"))

from lalib import prompt_evidence as PE  # noqa: E402


def test_opengenes_decrease_option_b():
    text = (
        "Animal model evidence:\nNo model organism data available for this gene.\n"
        "Functional annotations: cytoskeletal motor activity, microtubule binding, apoptotic process.\n"
        "In muscle of male humans between ages 21 and 67 years, Gene X mRNA expression:\n"
        "Options: A. increases with age B. decreases with age"
    )
    out = PE.run("opengenes-masked", text)
    assert out["status"] == "ok"
    assert out["suggested_letter"] in ("A", "B", None)
    # motor/apoptosis cues lean decrease → B when A=increase B=decrease
    assert out["option_map"]["decrease_letter"] == "B"


def test_synergyage_additive():
    text = (
        "Which intervention results in longer mean lifespan?\n"
        "Options: A. Pink1(KO) B. Pink1(KO);Snca(KO)\n"
        "Known individual effects:\n"
        "Pink1(KO): +8.6% lifespan (592 days vs 545 days WT) | "
        "Snca(KO): -39.2%, -39.2%. Mean: -39.2% (226 days vs 372 days WT)"
    )
    out = PE.run("synergyage-pairwise", text)
    assert out["status"] == "ok"
    assert out["suggested_letter"] == "A"
    assert out["confidence"] in ("high", "medium")


def test_cli_prompt_evidence():
    text = (
        "Options: A. Foo(KO) B. Bar(KO)\n"
        "Foo(KO): +10% lifespan | Bar(KO): -5% lifespan"
    )
    p = subprocess.run(
        [
            sys.executable,
            str(ROOT / "scripts" / "la.py"),
            "prompt-evidence",
            "--kind",
            "synergyage-pairwise",
            "--text",
            text,
        ],
        capture_output=True,
        text=True,
    )
    assert p.returncode == 0, p.stderr
    out = json.loads(p.stdout)
    assert out["suggested_letter"] == "A"
