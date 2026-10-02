"""Pea chromosomes reach ~650 Mb, beyond the 2^29 (536,870,912 bp) limit of
BAI indexes: every BAM index asmqc writes or lets CRAQ write must be CSI."""

import os
import shutil
import subprocess
from pathlib import Path

import pytest

from asmqc import mapping

SHIM = Path(__file__).parents[2] / "workflow" / "bin" / "craq_shim"
pytestmark = pytest.mark.skipif(shutil.which("samtools") is None, reason="no samtools")


@pytest.fixture
def big_bam(tmp_path) -> Path:
    sam = ("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr5\tLN:700000000\n"
           "r1\t0\tchr5\t600000000\t60\t10M\t*\t0\t0\tACGTACGTAC\tFFFFFFFFFF\n")
    bam = tmp_path / "big.bam"
    subprocess.run(["samtools", "view", "-b", "-o", str(bam), "-"], input=sam.encode(),
                   check=True)
    return bam


def test_bai_cannot_hold_pea_chromosomes(big_bam):
    r = subprocess.run(["samtools", "index", str(big_bam)], capture_output=True, text=True,
                       check=False)
    assert r.returncode != 0 and "csi" in r.stderr.lower()


def test_mapping_writes_csi(big_bam):
    mapping.index_bam(big_bam)
    assert Path(f"{big_bam}.csi").exists() and not Path(f"{big_bam}.bai").exists()
    out = subprocess.run(["samtools", "view", "-c", str(big_bam), "chr5:599999990-600000020"],
                         capture_output=True, text=True, check=True).stdout
    assert out.strip() == "1"


def test_craq_shim_turns_index_into_csi(big_bam):
    env = os.environ | {"PATH": f"{SHIM}:{os.environ['PATH']}"}
    subprocess.run(["samtools", "index", str(big_bam)], env=env, check=True)
    assert Path(f"{big_bam}.csi").exists()
    v = subprocess.run(["samtools", "--version"], env=env, capture_output=True, text=True,
                       check=True).stdout
    assert v.startswith("samtools")  # other subcommands pass through
